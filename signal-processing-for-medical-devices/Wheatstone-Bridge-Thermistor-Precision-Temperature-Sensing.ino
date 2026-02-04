/**************************************************************************

Wheatstone-Bridge-Thermistor-Precision-Temperature-Sensing

/**************************************************************************/

// Required libraries
#include <Wire.h>
#include <SPI.h>
#include <SD.h>

// Input parameters
float R1 = 10000 ;
float R2 = 10000 ;
float R3 = 10000;
float inVoltage = 3.3 ;
float sampRate = 1 ;  

// DAQ variables 
int V1;
int V2;

// Temp calculation variables
float Vdiff; // bridge node difference 
float resistance; // thermistor resistance 
float logR;

// Steinhart–Hart coefficients (from curve fitting)
float SHa = 1.13E-03;
float SHb = 2.35E-04;
float SHc = 8.54E-08;
float temp;


// Moving average variables
int i;
const int buffsize = 10; // window length of my moving average 
float tempBuff[buffsize]; // array tempBuff that can hold buffsize number of floats 
float sumTemp = 0; // running sum of the values currently in the buffer 
float avgTemp; // current moving average 
int bufIndex = 0; 


// Save data to SD card
const int chipSelect = 4;
File tempRecord;

// Initialize data acqusition and calculations
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

   tempRecord = SD.open("tempLog.txt", FILE_WRITE);
   if (tempRecord){
   tempRecord.println("Voltage (1)  Voltage (2)  Voltage (diff)  Resistance  Temp (C)  Average temp (C)");
   tempRecord.close();
   }
   else {
     Serial.println("error opening tempRecord.txt");
     return;
   }

  pinMode(A1, INPUT);
  pinMode(A2, INPUT);

}

void loop() {

  // Read voltage from Wheatstone bridge
  V1 = analogRead(A1);  
  V2 = analogRead(A2);

  // Calculate difference between nodes (equivalent to voltage across bridge)
  Vdiff = ((V1 - V2)*(3.3/1023)); // ADC  

  // Calculate resistance from voltage
  resistance = ((R2*inVoltage - (R1+R2) * Vdiff)/ ((R1*inVoltage) + (R1+R2) * Vdiff)) * R3; 

  // Calculate temperature from resistance
  logR = log(resistance);
  temp = 1/ (SHa + (SHb * logR) + (SHc * powf(logR,3)));
  temp = temp - 273.15;

  // Moving average Circular buffer 
  float oldest = 0.0;
  oldest = tempBuff[bufIndex];
  sumTemp = sumTemp - oldest;
  tempBuff[bufIndex] = temp;
  sumTemp = sumTemp + temp;

  // Advance the circular index - wrap at buffsize so it goes back to 0
  bufIndex = bufIndex + 1;
  if (bufIndex == buffsize) bufIndex = 0;

  // Compute the average over the current window (grows until full)
  avgTemp = sumTemp/buffsize;


  // Output (Serial.print prints to the serial monitor; tempRecord.print prints to the SD card)
  Serial.print(V1);
  Serial.print("\t");
  Serial.print(V2);
  Serial.print("\t");
  Serial.print(Vdiff);
  Serial.print("\t");
  Serial.print(resistance);
  Serial.print("\t");
  Serial.print(temp);
  Serial.print("\t");
  Serial.println(avgTemp);

  // Uncomment to use SD card
  tempRecord = SD.open("tempLog.txt", FILE_WRITE);
  
   if (tempRecord) {
     tempRecord.print(V1);
     tempRecord.print("\t");
     tempRecord.print(V2);
     tempRecord.print("\t");
     tempRecord.print(Vdiff);
     tempRecord.print("\t");
     tempRecord.print(resistance);
     tempRecord.print("\t");
     tempRecord.print(temp);
     tempRecord.print("\t");
     tempRecord.println(avgTemp);
     tempRecord.close();
    }
   else {
     //Serial.println("error opening tempLog.txt");
     return;
   }
  
  // Implement sampling rate
  delay(1000);  

}



