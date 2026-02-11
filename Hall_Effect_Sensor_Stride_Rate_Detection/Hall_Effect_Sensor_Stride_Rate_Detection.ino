/**************************************************************************
Lab 5 
Hall Effect Sensor - Stride Rate Detection
Signal Processing for Medical Devices 2025

This project implements real-time stride-rate estimation using a Hall effect sensor and a magnet to detect periodic knee motion (prosthesis-inspired gait sensing concept). The pipeline converts raw ADC readings into magnetic field strength (Gauss), detects stride events in the time domain, smooths cadence estimates using a moving average, and independently estimates cadence via FFT.

## What this measures
- A Hall effect sensor mounted near a rotating joint measures magnetic field changes as a magnet passes close to the sensor.
- The resulting waveform behaves like a switch-like square wave: one saturated field value when the magnet is away and the opposite polarity when the magnet is within trigger range.
- We extract cadence as:
  1) event-detected stride intervals (time-domain),
  2) smoothed event cadence (moving average),
  3) dominant periodicity via FFT (frequency-domain). :contentReference[oaicite:3]{index=3}

---------------------- APPLICATIONS ----------------------
• Wearable gait monitoring
• Prosthetic/rehabilitation sensing
• Embedded physiological signal processing
• Real-time cadence estimation from magnetic sensors

---

## Hardware Components
- Arduino (10-bit ADC)
- Analog Hall effect sensor (output voltage centered around Vcc/2)
- Small permanent magnet

### Wiring (typical)
- Hall sensor Vcc → 3.3V (or 5V if your sensor supports it; update `inVoltage` / `zeroFieldVoltage`)
- Hall sensor GND → GND
- Hall sensor OUT → A1
- SD module CS → D4 (matches `chipSelect = 4`)
- SD module MOSI/MISO/SCK → Arduino SPI pins


## Signal pipeline
### 1) ADC → Voltage
We convert raw 10-bit ADC reading to volts:
Vout = rawValue / (2^10 − 1) * inVoltage

### 2) Voltage → Magnetic Field (Gauss)
magneticField = (Vout − zeroFieldVoltage) / (sensitivity in V/G)

### 3) Event Detection (time-domain)
A stride event is detected when:
- magneticField < MAG_FIELD_THRESHOLD
- derivative < DERIV_THRESHOLD
- AND debounce time has passed (prevents double-triggers)

This outputs stride interval (s) and cadence (steps/min) = 60 / stride_interval.

### 4) Moving Average (cadence smoothing)
A circular buffer (N=5) smooths cadence to reduce jitter and step-to-step variability, producing a more clinically meaningful cadence trend. 

### 5) FFT cadence estimate (frequency-domain)
A 256-sample window is mean-centered, windowed (Hamming), FFT’d, and the dominant frequency between 0.3–5 Hz is selected and converted to steps/min.

With sampRate=50 Hz and N=256:
- frequency resolution df = Fs/N = 50/256 = 0.1953 Hz
- this corresponds to 11.7 steps/min resolution (0.1953*60)

 **************************************************************************/

// Required libraries
#include <Wire.h>
#include <SPI.h>
#include <SD.h>
#include <arduinoFFT.h>

// Input parameters
float inVoltage = 3.3;         // Reference voltage (V)
float sampRate = 50;           // Sampling rate (Hz) for Hall sensor
float sensitivity = 5.0;       // Hall sensor sensitivity (mV/G) 
float zeroFieldVoltage = 1.65; // Output voltage at zero magnetic field (V) = Vcc/2

// Event detection thresholds
const float MAG_FIELD_THRESHOLD = 0.0;      // Magnetic field threshold for event (units = Gauss)
const float DERIV_THRESHOLD = -300;         // Derivative threshold for event detection
const unsigned long DEBOUNCE_TIME = 300;    // Minimum time between events (ms) to avoid noise due to repeated triggers - set to 300ms as this is the lower bound of time between strides 
const unsigned long EVENT_TIMEOUT = 10000;  // Reset stride detection after 5 seconds of inactivity

// DAQ variables
int rawValue;                 // Raw ADC reading
float Vout;                   // Output voltage from Hall sensor (V)
float magneticField;          // Magnetic field strength (Gauss)
float prev_magneticField = 0; // Previous magnetic field value
float deriv_magneticField;    // Derivative of magnetic field
int gaitState = 0;            // When gaitState = 1, this indicates that an overlap event has been detected and the magnet is within the trigger range relative to the sensor. 

// Event detection & stride rate variables
unsigned long event_time = 0;      // Time of current event 
unsigned long prev_event_time = 0; // Time of previous event 
const int max_strides = 10;        // Size of the history buffer for stride intervals
float stride_times[max_strides];   // Array to store stride intervals
int stride_count = 0;              // Number of recorded stride events
float current_stride_rate = 0;     // Current calculated stride rate (steps/min)

// FFT variables for frequency analysis
const int FFT_SAMPLES = 256;     
double fftReal[FFT_SAMPLES];         // Real array for FFT input/output 
double fftImag[FFT_SAMPLES];         // Imaginary array for FFT input/output 
ArduinoFFT<double> gaitFFT(fftReal, fftImag, FFT_SAMPLES, sampRate); 
int fftIndex = 0;                    // Index into the FFT buffer
double fft_stride_rate = 0.0;        // Stride rate calculated from FFT in steps/min 
void performGaitFFT();               // Function declaration so it can be called in loop (). Full function is at the end of code

// Moving average variables for stride times and stride rate smoothing 
int i = 0;
const int STRIDE_AVG_BUFFSIZE = 5;  // Buffer for smoothing stride intervals
float sum_stride_time;              // Sum of stride times 
float avg_stride_time;              // Average of stride times 
float stride_interval_buff[STRIDE_AVG_BUFFSIZE];
float stride_interval_sum = 0;
int stride_interval_index = 0;
float smoothed_stride_rate = 0;  // Smoothed stride rate output

// Save data to SD card 
const int chipSelect = 4;
File hallRecord;

// Initialize data acquisition and calculations
void setup() {
  Serial.begin(9600);
  while (!Serial);
  

  //Uncomment to use SD card
  SD.begin(chipSelect);
  if (!SD.begin(chipSelect)) {
    Serial.println("SD not found ");
    return;
  }
  Serial.println("SD found");
  
  hallRecord = SD.open("hallLog.txt", FILE_WRITE);
  if (hallRecord) {
    hallRecord.println("Raw Value  Voltage (V)  Magnetic Field (G)  Average Field (G)");
    hallRecord.close();
  }
  else {
    Serial.println("error opening hallLog.txt");
    return;
  }
  
  pinMode(A1, INPUT);
}

void loop() {

  // Voltage -> Magnetic field 
  // Read voltage from Hall effect sensor on pin A1
  rawValue = analogRead(A1);
  
  // Convert ADC reading to voltage
  Vout = rawValue / (pow(2, 10) - 1) * inVoltage;
  
  // Calculate magnetic field from voltage: Field (Gauss) = (Vout - ZeroFieldVoltage) / (Sensitivity in V/G)
  prev_magneticField = magneticField; 
  magneticField = (Vout - zeroFieldVoltage) / (sensitivity / 1000.0); 

  //****** 1. Event Detection *******
  // A sensor/magnet overlap event is defined as 1) negative magnetic field 2) derivative < 300 (B = -328)
  deriv_magneticField = magneticField - prev_magneticField;

  // Detect stride event if conditions are obeyed 
  if (magneticField < MAG_FIELD_THRESHOLD && deriv_magneticField < DERIV_THRESHOLD) {
    unsigned long current_time = millis(); 
    
    // Check debounce time to avoid multiple triggers from same magnet pass
    if (current_time - prev_event_time > DEBOUNCE_TIME) {
      gaitState = 1;   // magnet close so joint overlapping = event is triggered
      
      // Calculate stride time interval
      if (prev_event_time > 0) {
        float stride_interval = (current_time - prev_event_time) / 1000.0; // Convert to seconds
        Serial.print("stride_interval:");
        Serial.print(stride_interval);
        Serial.print("\t");
        
        // ****** 2.1. Moving average of stride intervals *******
        if (stride_count < max_strides) {
          stride_times[stride_count] = stride_interval;
          stride_count++;
        } else {
          
          for (int j = 0; j < max_strides - 1; j++) {
            stride_times[j] = stride_times[j + 1];
          }
          stride_times[max_strides - 1] = stride_interval;
        }

        sum_stride_time = 0;
        for (int j = 0; j < max_strides; j++) {
           sum_stride_time += stride_times[j];
        }
        avg_stride_time = sum_stride_time / stride_count; 
      
        // ****** 2.2. Moving average of stride rates *******
        float stride_rate_from_interval = 60.0 / stride_interval; // Compute instantaneous stride rate
        stride_interval_sum = stride_interval_sum - stride_interval_buff[stride_interval_index];
        stride_interval_buff[stride_interval_index] = stride_rate_from_interval;
        stride_interval_sum = stride_interval_sum + stride_interval_buff[stride_interval_index];
        stride_interval_index++;
        stride_interval_index %= STRIDE_AVG_BUFFSIZE;

        smoothed_stride_rate = stride_interval_sum / STRIDE_AVG_BUFFSIZE;

      }
      
      prev_event_time = current_time;
      event_time = current_time;
    }
  } 

    
  else {
  // Event Timeout  - this resets stride detection after 5 seconds of inactivity e.g if person is standing still 
    if (millis() - prev_event_time > EVENT_TIMEOUT && prev_event_time > 0) {
    gaitState = 0;
    stride_count = 0;
    } 
    else if (prev_event_time == 0) {
    gaitState = 0;
    }
  }

  // Store filtered signal for FFT analysis
  if (fftIndex < FFT_SAMPLES) {
    fftReal[fftIndex] = (double)magneticField; // REAL 
    fftImag[fftIndex] = 0.0; // IMAGINARY 
    fftIndex++;
  }

  // Run the FFT function below once buffer is full
    if (fftIndex >= FFT_SAMPLES) {
      performGaitFFT();  
      fftIndex = 0;       // reset for next window
    }


  
  // Output to Serial Monitor
  Serial.print(magneticField, 2);
  Serial.print("\t");
  Serial.print(avg_stride_time);
  Serial.print("\t");
  Serial.print(smoothed_stride_rate, 1); // **** removed the stride rate moving average 
  Serial.print(" steps/min");
  Serial.print("\t");
  Serial.print("FFTGaitRate: ");
  Serial.print(fft_stride_rate, 1);
  Serial.println(" steps/min");


  hallRecord = SD.open("hallLog.txt", FILE_WRITE);
  if (hallRecord) {

    hallRecord.print(magneticField, 2);
    hallRecord.print("\t");

    hallRecord.print(avg_stride_time, 3);
    hallRecord.print("\t");

    hallRecord.print(smoothed_stride_rate, 1);
    hallRecord.print("\t");

    hallRecord.print(fft_stride_rate, 1);

    hallRecord.println();   
    hallRecord.close();
  }
  delay(1000 / sampRate);


}


//****** 3. FFT *******
void performGaitFFT() {
  
  // Remove DC offset
  double meanB = 0.0;
  for (int i = 0; i < FFT_SAMPLES; i++) {
    meanB += fftReal[i];
  }
  meanB /= (double)FFT_SAMPLES;
  for (int i = 0; i < FFT_SAMPLES; i++) {
    fftReal[i] -= meanB;
  }

  // Apply window to reduce spectral leakage
  gaitFFT.windowing(FFTWindow::Hamming, FFTDirection::Forward);
  // Compute FFT
  gaitFFT.compute(FFTDirection::Forward); 
  // Convert complex spectrum to magnitude spectrum in fftReal[]
  gaitFFT.complexToMagnitude();

  // Find dominant gait frequency in a given range 0.3–5 Hz - Context: these are upper and lower bounds for stride rate in steps/second
  double minFreq = 0.3; // walking lower bound
  double maxFreq = 5.0; // upper bound (fast run)

  // Convert to FFT bins 
  int minBin = (int)(minFreq * FFT_SAMPLES / sampRate);  
  int maxBin = (int)(maxFreq * FFT_SAMPLES / sampRate);

  // Supplementary conditions  
  if (minBin < 1) minBin = 1;                           // ignore DC
  if (maxBin > FFT_SAMPLES / 2) maxBin = FFT_SAMPLES/2; // limit to Nyquist  

  double peakMag = 0.0;
  int peakBin = 0;
  for (int k = minBin; k <= maxBin; k++) {
    if (fftReal[k] > peakMag) {
      peakMag = fftReal[k];
      peakBin = k;
    }
  }

  // Convert peak bin to frequency value in Hz (steps/second)
  double fpeak = peakBin * sampRate / FFT_SAMPLES;              

  // Convert to steps/min 
  double stepsPerMin = fpeak * 60.0;
  fft_stride_rate = stepsPerMin;


}




