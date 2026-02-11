 /*
  Photoplethysmography_FrequencyDomain_HeartRate_FFT.ino

  Summary
  -------
  Estimates heart rate (BPM) from MAX30102 photoplethysmography using FFT-based
  frequency-domain peak detection on the IR channel.

  Pipeline
  --------
  1) Acquire N samples (IR + Red) from MAX30102
  2) Remove DC offset from IR (zero-mean)
  3) Apply Hamming window to reduce spectral leakage
  4) Compute FFT and magnitude spectrum
  5) Search for dominant peak within a physiological HR band (e.g., ~0.7–3.3 Hz)
  6) Convert peak frequency to BPM and smooth with moving-average
  7) Optional: compare against Maxim library HR/SpO2 outputs for validation

  Notes
  -----
  - samplingFrequency MUST match the true sample rate used for the FFT window.
  - Searching only within a physiological band reduces false peaks from drift/respiration.
  - Peak loc is hitting 10 as a minimum - this is quite high (60bpm as minimum bpm) 
   However I determined that the removal of the DC offset in the code was not sufficient without external filtering
   there could be a peak that lies after the actual 0 peak that is higher than the 0 peak - hence why i was getting previously low values 
*/


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

// sensor input/output variables
int i;
uint32_t irBuffer[256]; 
uint32_t redBuffer[256];  
byte readLED = 13; 
byte ledBrightness = 60;  //Options: 0=Off to 255=50mA
byte sampleAverage = 4;  //Options: 1, 2, 4, 8, 16, 32
byte ledMode = 2;       //Options: 1 = Red only, 2 = Red + IR, 3 = Red + IR + Green
int sampleRate = 200;  //Options: 50, 100, 200, 400, 800, 1000, 1600, 3200
int pulseWidth = 411; //Options: 69, 118, 215, 411
int adcRange = 4096; //Options: 2048, 4096, 8192, 16384

// variables for maxim SpO2/HR calculations
int32_t bufferLength = 256; 
int32_t spo2;
int8_t validSPO2; 
int32_t heartRate; 
int8_t validHeartRate; 

// Variables to calculate HR
const uint16_t samples = 256; 
double vReal[samples];
double vImag[samples];
float samplingFrequency = 25; 
float irnoDC;

ArduinoFFT<double> FFT = ArduinoFFT<double>(vReal, vImag, samples, samplingFrequency);

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
float irDC;
float myavgHR;
int j;

void setup() {

  Serial.begin(115200);
  Serial.println("Initializing...");

//  pinMode(pulseLED, OUTPUT);
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
    {
      PPGSensor.check(); 
    }
    redBuffer[i] = PPGSensor.getRed();
    irBuffer[i] = PPGSensor.getIR();
    PPGSensor.nextSample(); 
  }

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

      PPGSensor.nextSample(); 
    }


  
  // Calculate HR and SpO2 using existing library
  maxim_heart_rate_and_oxygen_saturation(irBuffer, bufferLength, redBuffer, &spo2, &validSPO2, &heartRate, &validHeartRate);
  
  // Remove DC component 
  float irnoDC[samples];
  float meanIR = 0.0;
  for (int i = 0; i < samples; i++){
    meanIR += irBuffer[i];
  }
  meanIR /= samples;
  for (int i = 0; i < samples; i++) {  
  irnoDC[i] = (float)irBuffer[i] - meanIR;   
  }


  //double ratio = twoPi * signalFrequency / samplingFrequency; Fraction of a complete cycle stored at each sample (in radians)
  for (uint16_t i = 0; i < samples; i++)
  {
    vReal[i] = irnoDC[i]; /* Build data with positive and negative values
    //vReal[i] = uint8_t((amplitude * (sin(i * ratio) + 1.0)) / 2.0);/* Build data displaced on the Y axis to include only positive values*/
    vImag[i] = 0.0; //Imaginary part must be zeroed in case of looping to avoid wrong calculations and overflows
  }
  FFT.windowing(FFTWindow::Hamming, FFTDirection::Forward);	  /* Weigh data */
  FFT.compute(FFTDirection::Forward);                        /* Compute FFT */
  FFT.complexToMagnitude();                              /* Compute magnitudes */
  //PrintVector(vReal, samples>>1, SCL_PLOT); // Will print magnitude of all peaks
  

  // Find index of maximum peak
  // Don’t want the algorithm to pick breathing (~0.2–0.4 Hz) or motion artifacts & Heart rate lies between 0.7–3.3 Hz (approx 42–200 bpm). 
  // So only search only those FFT bins. 
  int kMIN = (60.0 / 60.0) * samples / samplingFrequency; 
  int kMAX = (200.0 / 60.0) * samples / samplingFrequency;

  double peak = -1.0; 
  peakloc = 0;
  for ( int i = kMIN; i <= kMAX; i++) {
    if (vReal[i] > peak) {
      peak = vReal[i];
      peakloc = i;
    }
  }

  // Convert index of max peak to frequency 
  float fpeak = 0.0;
  fpeak = peakloc * samplingFrequency / samples;

  // Convert frequency to HR 
  myHR = 60.0 * fpeak;

  // Calculate moving average HR 
  HRBuff[j] = myHR;
  j = (j + 1) % buffsize;
  sumHR = 0;
  for (int j = 0; j < buffsize; j++) {
    sumHR += HRBuff[j];
  }
  myavgHR = sumHR / buffsize;



    
    // Serial.print(spo2, DEC);
    // Serial.print("\t");
    // Serial.print(mySpO2, 0);
    // Serial.print("\t");
    Serial.print(heartRate,DEC);
    Serial.print("\t");
    Serial.print(myHR,0);
    Serial.print("\t");
    Serial.print(myavgHR,0);
    Serial.println("\t");
   // Serial.print(peakloc, 0);
    //Serial.print("\t");
    //Serial.print(kMIN, 0);
    //Serial.print("\t");
    // Serial.print(kMAX, 0);
    //Serial.println("\t");
    
  // pulseRecord = SD.open("HRLog.txt", FILE_WRITE);
  
  // if (pulseRecord) {
  //   pulseRecord.print(heartRate, DEC);
  //   pulseRecord.print("\t");
  //   pulseRecord.print(myavgHR);
  //   pulseRecord.print("\t");
  //   pulseRecord.close();
  // }
  // else {
   // Serial.println("error opening HRLog.txt");
  //   return;
  // }
  
  }
    
}


