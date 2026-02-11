/**************************************************************************

Photoplethysmography_Autocorrelation_HeartRate.ino

  Summary
  -------
  Estimates heart rate from reflective photoplethysmography (PPG) acquired with a MAX30102.
  Uses autocorrelation on a DC-removed IR waveform to identify the dominant cardiac period
  (lag), then converts lag (samples) to BPM. Also computes a moving-average smoothed HR.

  Signal Processing Pipeline
  --------------------------
  1) Acquire a fixed window of PPG samples (Red + IR) from MAX30102
  2) Remove DC offset from IR (zero-mean normalization)
  3) Autocorrelation over physiologic lag bounds (e.g., 40–150 BPM mapped to sample lags)
  4) Select best lag (maximum autocorrelation peak) → HR = sampleRate * 60 / lag
  5) Smooth HR with moving average for stability
  6) Optional: compare against Maxim library HR/SpO2 outputs; log to Serial/SD as CSV

  Key Parameters
  --------------
  - sampleRate:        200 Hz (sensor setting)
  - window length:     256 samples (power-of-2; supports FFT in other sketches)
  - HR bounds:         lag search constrained to physiologic range to reject artifacts


/**************************************************************************/

#include <Wire.h>
#include "MAX30105.h"
#include "heartRate.h"
#include "spo2_algorithm.h"
#include "arduinoFFT.h"
#include <SPI.h>
#include <SD.h>

MAX30105 PPGSensor;

#define I2C_SPEED_FAST 400000
#define MAX_BRIGHTNESS 255

const int chipSelect = 4;
File pulseRecord;

// Sensor input/output variables
int i;
uint32_t irBuffer[256]; //infrared LED sensor data
uint32_t redBuffer[256];  //red LED sensor data
byte readLED = 13; //Blinks with each data read
byte ledBrightness = 60; //Options: 0=Off to 255=50mA
byte sampleAverage = 4; //Options: 1, 2, 4, 8, 16, 32
byte ledMode = 2; //Options: 1 = Red only, 2 = Red + IR, 3 = Red + IR + Green
int sampleRate = 200; //Options: 50, 100, 200, 400, 800, 1000, 1600, 3200
int pulseWidth = 411; //Options: 69, 118, 215, 411
int adcRange = 4096; //Options: 2048, 4096, 8192, 16384

// Variables for pre-existing maxim SpO2/HR calculations
int32_t bufferLength = 256; //changed from 256 for moving average to compute faster
int32_t spo2; //SPO2 value
int8_t validSPO2; //indicator to show if the SPO2 calculation is valid
int32_t heartRate; //heart rate value
int8_t validHeartRate; //indicator to show if the heart rate calculation is valid

// Variables to calculate HR
const uint16_t samples = 256; 
float sum;
float irDC;
#define SCL_INDEX 0x00
#define SCL_TIME 0x01
#define SCL_FREQUENCY 0x02
#define SCL_PLOT 0x03
double peak;
float peakloc;
float currIR;
float diffIR[256];
float myHR;
const int buffsize = 10;
float HRBuff[buffsize];
float sumHR;
float myavgHR;
int j;
float lag;

void setup() {

  Serial.begin(115200);
  Serial.println("Initializing...");

 // pinMode(pulseLED, OUTPUT);
  pinMode(readLED, OUTPUT);

  SD.begin(chipSelect);
  if (!SD.begin(chipSelect)) {
    Serial.println("SD not found ");
    return;
  }
  Serial.println("SD found");
  File pulseRecord = SD.open("HRLog.txt", FILE_WRITE);
  if (pulseRecord){
  pulseRecord.println("SpO2 (library)  SpO2 (my calcs)  PR (library)  PR (my calcs)  HRV");
  pulseRecord.close();
  }
  else {
    Serial.println("error opening HRLog.txt");
    return;
  }

  // Initialize sensor
  if (!PPGSensor.begin(Wire, I2C_SPEED_FAST)) {
    Serial.println("MAX30102 was not found. Please check wiring/power. ");
    while (1);
  }

  PPGSensor.setup(ledBrightness, sampleAverage, ledMode, sampleRate, pulseWidth, adcRange); //Configure sensor with these settings
  PPGSensor.setPulseAmplitudeRed(0x0A); //Turn Red LED to low to indicate sensor is running
  PPGSensor.setPulseAmplitudeGreen(0); //Turn off Green LED
}

void loop() {
  
  // Read the first 100 samples, and determine the signal range
  for (i = 0 ; i < bufferLength ; i++)
  {
    while (PPGSensor.available() == false) 
      PPGSensor.check(); 

    redBuffer[i] = PPGSensor.getRed();
    irBuffer[i] = PPGSensor.getIR();
    PPGSensor.nextSample(); 
  }


  // Continuously taking samples from MAX30102. Heart rate and SpO2 are calculated every 1 second
  while (1)
  {
    // Dumping the first 25 sets of samples in the memory and shift the last 75 sets of samples to the top
    for (int i = 64; i < 256; i++)
    {
      redBuffer[i - 64] = redBuffer[i];
      irBuffer[i - 64] = irBuffer[i];
      diffIR[i - 64] = diffIR[i];
    }

    // Take 25 sets of samples before calculating the heart rate.
    for (int i = 192; i < 256; i++)
    {
      while (PPGSensor.available() == false) 
      {
        PPGSensor.check();
      } 

      digitalWrite(readLED, !digitalRead(readLED));

      redBuffer[i] = PPGSensor.getRed();
      irBuffer[i] = PPGSensor.getIR();
       float prevIR = (i > 0) ? (float)irBuffer[i - 1] : (float)irBuffer[i];;
      currIR = (float)irBuffer[i]; 
      diffIR[i] = currIR - prevIR;

      PPGSensor.nextSample(); 
    }
  
    //  --- Library HR/SpO2 (validation) ---
    maxim_heart_rate_and_oxygen_saturation(irBuffer, bufferLength, redBuffer, &spo2, &validSPO2, &heartRate, &validHeartRate);
  
    // --- DC removal ---
    float zeromeanIR[256];             
    float meanIR = 0.0f;
    for (int i = 0; i < samples; i++) {
      meanIR += (float)irBuffer[i];       
      }
      meanIR /= (float)samples;          
    for (int i = 0; i < samples; i++) {
      zeromeanIR[i] = (float)irBuffer[i] - meanIR;
      }

    //  --- Autocorrelation HR estimation ---

    // Converting HR bounds to lag bounds in samples: define a search band of lags that lie between an acceptable heart rate range (150 - 40 bpm) 
    int MINlag = sampleRate * 60 / 150; // (80 samples)
    int MAXlag = sampleRate * 60 / 40; // (300 samples)
    if (MAXlag > (samples - 1)) MAXlag = samples - 1; 
    
    // Storing the best lag with the largest autocorrelation in "peak" and track its location in "peakloc" 
    float peak  = -1e30f;  
    float peakloc  = -1;  

    for (int lag = MINlag; lag <= MAXlag; lag++) {
      float sum = 0.0;
      for (int k = 0; k < samples - lag; k++) {
      sum += zeromeanIR[k] * zeromeanIR[k + lag];
    }

    // Find index of maximum value 
    if (sum > peak) {
    peak = sum;
    peakloc = lag;
    }
    }
  
    // Clamp to window so the inner loop has work to do
    if (MAXlag >= samples) MAXlag = samples - 1; 
    if (MINlag < 1)        MINlag = 1;

    // Convert index of max value to HR 
    myHR = sampleRate * 60 / peakloc;

    // Calculate moving average HR 
    HRBuff[j] = myHR;
    j = (j + 1) % buffsize;
    sumHR = 0;
    for (int n = 0; n < buffsize; n++) {
      sumHR += HRBuff[n];
    }
    myavgHR = sumHR / buffsize;
 

    // Serial.print(spo2, DEC);
    // Serial.print("\t");
    // Serial.print(mySpO2, 0);
    // Serial.print("\t");
    Serial.print(heartRate, DEC);
    Serial.print("\t");
    Serial.print(myHR, 0);
    Serial.print("\t");
    Serial.print(myavgHR, 0);
    Serial.print("\t");
    Serial.print(peakloc, 0);
    Serial.println("\t");
  }
}

  
