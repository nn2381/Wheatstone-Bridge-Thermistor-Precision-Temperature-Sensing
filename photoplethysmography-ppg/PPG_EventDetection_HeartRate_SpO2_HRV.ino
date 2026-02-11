
/*
  Photoplethysmography_EventDetection_HeartRate_SpO2_HRV.ino

  Summary
  -------
  Reads reflective photoplethysmography (PPG) from a MAX30102 (Red + IR) and estimates:
    - Heart Rate (BPM): time-domain event detection on IR using derivative sign change
    - HRV: variability from inter-beat intervals (IBI) using running variance / SDNN-style metric
    - SpO2: manual ratio-of-ratios using AC/DC components of Red and IR

  Algorithm (Time Domain)
  -----------------------
  1) Sliding window acquisition (bufferLength samples) from MAX30102
  2) Derivative of IR: diffIR[n] = IR[n] - IR[n-1]
  3) Peak detection when slope transitions from + to − (crossedPosToNeg),
     gated by amplitude threshold (IRthresh) and refractory period
  4) IBI computed from sample index differences; HR = 60 / IBI_sec
  5) Moving-average smoothing of HR
  6) SpO2 computed from AC/DC ratio-of-ratios:
       R = (AC_red/DC_red) / (AC_ir/DC_ir),  SpO2 ≈ 110 − 25R
  7) Optional: compare against Maxim library HR/SpO2 for validation

  Notes
  -----
  - This sketch is designed for clarity and validation rather than optimized low-power deployment.
*/


#include <Wire.h>
#include "MAX30105.h" 
#include "heartRate.h"
#include "spo2_algorithm.h"

MAX30105 PPGSensor;

#define I2C_SPEED_FAST 400000
#define MAX_BRIGHTNESS 255


// Sensor input/output variables
int i;
uint32_t irBuffer[200]; 
uint32_t redBuffer[200];  
byte readLED = 13; 
byte ledBrightness = 60; //Options: 0=Off to 255=50mA
byte sampleAverage = 4; //Options: 1, 2, 4, 8, 16, 32
byte ledMode = 2;      //Options: 1 = Red only, 2 = Red + IR, 3 = Red + IR + Green
int sampleRate = 200;
int pulseWidth = 411; //Options: 69, 118, 215, 411
int adcRange = 4096; //Options: 2048, 4096, 8192, 16384

// Variables for maxim SpO2/HR calculations
int32_t bufferLength = 200; 
int32_t spo2; 
int8_t validSPO2; 
int32_t heartRate; 
int8_t validHeartRate; 

// Variables to calculate SpO2
float maxRed;
float minRed;
float maxIR;
float minIR;
float sumRed;
float sumIR;
float redDC;
float irDC;
float redAC;
float irAC;
float R;
float mySpO2;

// Variables to calculate HR
float currIR;
float diffIR[200];
byte isPos = 0;
byte prevPos = 0;
float diffThresh = 8.0;
float IRthresh = 950.0 ;
int beatDetected;
float beatTime;
float prevBeatTime;
float timeBetweenBeats;
float myHR;
const int buffsize = 10;
float HRBuff[buffsize];
float sumHR;
float myavgHR;
int j;
const float fs = 25.0f;              // ACTUAL hidden sampling rate in the sensor 
static long globalIdx = 0;          // total samples since start 
static long prevBeatIdx = -100000; // sample index of last beat
float  IBI_sec = 0;               // interbeat interval in seconds 

// Variables for HRV calculation
float IBIavg = 0.0;
float IBI = 0.0;
float sumVar = 0.0;
float variance = 0.0;
float HRV = 0.0;
int IBICount = 0;

void setup() {

  Serial.begin(115200);
  Serial.println("Initializing...");

  pinMode(readLED, OUTPUT);

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
  
  // Creating a window of length 200 samples
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

  //Sliding the window by 50 samples
  while (1)
  {
    //dumping the first 50 sets of samples in the memory and shift the last 150 sets of samples to the top
    for (byte i = 50; i < 200; i++)
    {
      redBuffer[i - 50] = redBuffer[i];
      irBuffer[i - 50] = irBuffer[i];
      diffIR[i - 50] = diffIR[i];
    }

    //Adding the new 50 samples to the sliding window
    for (byte i = 150; i < 200; i++)
    {
      while (PPGSensor.available() == false) 
      {
        PPGSensor.check(); 
      }

      digitalWrite(readLED, !digitalRead(readLED)); 

      redBuffer[i] = PPGSensor.getRed();
      irBuffer[i] = PPGSensor.getIR();

      // Calculate the derivative for event detection 
      float prevIR = (float)irBuffer[i - 1];
      currIR       = (float)irBuffer[i];
      diffIR[i]    = currIR - prevIR;

      // Zero-crossing detection
      prevPos = isPos;
      if      (diffIR[i] >  diffThresh) isPos = 1;   // slope is definitely rising
      else if (diffIR[i] < -diffThresh) isPos = 0;   // slope is definitely falling

      // Calculate HR manually using event detection //HR detection uses the IR component 
      bool crossedPosToNeg = (prevPos == 1) && (isPos == 0);  // Detect peak when slope flips from + to −
      bool highAmplitude   = (currIR > IRthresh);            // Check that the amplitude high enough
      int refractorySamples = (int)(0.4f * fs);             //Check that the refractory period in samples is reasonable -  0.30 s * fs samples must pass since last beat
      bool refractoryOK     = ((long)globalIdx - prevBeatIdx) > refractorySamples;

      if (crossedPosToNeg && highAmplitude && refractoryOK) {
        beatDetected = 1;

        long beatIdx = globalIdx;                        // sample index of this beat 

        if (prevBeatIdx > 0) {
          long  IBI_samples = beatIdx - prevBeatIdx;   // change in samples
          IBI_sec     = IBI_samples / fs;             // change in seconds

          if (IBI_sec > 0.30f && IBI_sec < 2.00f) {    
            myHR = 60.0f / IBI_sec;                    
          }
        }

        prevBeatIdx = beatIdx;

        // Calculate moving average HR 
        HRBuff[j] = myHR;
        j = (j + 1) % buffsize;
        sumHR = 0;
        for (int n = 0; n < buffsize; n++) {
          sumHR += HRBuff[n];
        }
        myavgHR = sumHR / buffsize;


      // Calculate HRV manually 
      IBICount++;                                               // Counting the number of valid IBI
      IBI = IBI_sec;                         
      float delta  = IBI - IBIavg;                             // Compute the deviation of this new IBI from the previous running mean
      IBIavg = ((IBIavg * (IBICount - 1)) + IBI) / IBICount;  // Update the running mean to include this new value
      float delta2 = IBI - IBIavg;                           // The new deviation of this IBI from the updated mean.
      sumVar += delta * delta2;                             // Update the running sum of squared deviations 
      variance = (sumVar / (IBICount - 1)); 
      HRV = sqrtf(variance);                 
      
      

        
        Serial.print("My HR = "); Serial.print(myHR, 1);
        Serial.print("\t");
        Serial.print("Average HR = "); Serial.print(myavgHR, 1);
        Serial.print("\t");
        Serial.print("My HRV = "); Serial.print(HRV, 1);
        Serial.println("\t");

 
      }

      PPGSensor.nextSample();   
      globalIdx++;              
    }
    // Calculate HR and SpO2 using existing library - For comparison to your calculated values
    maxim_heart_rate_and_oxygen_saturation(irBuffer, bufferLength, redBuffer, &spo2, &validSPO2, &heartRate, &validHeartRate);
      
    // Calculate SpO2 manually using equation from TI application report 
    float maxRed = 0.0f; 
    float minRed = 1000000.0f;
    float sumRed = 0.0f;

    float maxIR = 0.0f; 
    float minIR = 1000000.0f;
    float sumIR = 0.0f;

    for (int n = 0; n<bufferLength; n++) { 
        float red = redBuffer[n]; 
        float ir = irBuffer[n];
        if (red < minRed) minRed = red; 
        if (red > maxRed) maxRed = red;  
        if (ir < minIR)   minIR = ir; 
        if (ir > maxIR)   maxIR = ir; 
        sumRed += red; 
        sumIR += ir; 
      }

    // DC is estimated by mean & AC by peak-to-peak (max - min)
    float mySpO2 = 0;
    redDC = sumRed / bufferLength; 
    redAC = (maxRed - minRed)/2;

    irDC = sumIR / bufferLength; 
    irAC = (maxIR - minIR)/2; 

    R = (redAC / redDC) / ( irAC / irDC ); 
    mySpO2 = 110 - 25 * R;
    
    Serial.print("\t");
    Serial.print("My SpO2 = "); Serial.print(mySpO2, 1);
    Serial.println("\t");
  }
}



