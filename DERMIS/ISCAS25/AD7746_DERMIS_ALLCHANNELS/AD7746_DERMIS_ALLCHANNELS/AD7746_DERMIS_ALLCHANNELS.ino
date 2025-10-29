#include <Wire.h>
#include "AD7746.h"

// Initialize AD7746 object
AD7746 ad7746;


void setup() {
  // Initialize I2C communication
  Wire.begin();
  
  // Initialize serial communication
  Serial.begin(9600);
  
  // Initialize AD7746
  ad7746.initialize();
  
  // Check connection
  if (ad7746.testConnection()) {
    Serial.println("AD7746 connection successful");
  } else {
    Serial.println("AD7746 connection failed");
    while (1); // Loop forever if connection fails
  }

  // around +4.5pF parasitic offset capacitance. calibrated with CAPDACA=4pF and offset register (0.5pF)


  // Reset AD7746 device and configure registers
  ad7746.reset();
  //ad7746.writeCapSetupRegister(0xA0);
  //ad7746.writeCapSetupRegister(0xC0); // CAPEN=1, CN2=1, DIFF=0, rest=0
  ad7746.writeCapSetupRegister(0x80); // CAPEN=1, CN2=0, DIFF=0, rest=0
  //ad7746.writeExcSetupRegister(0x0B); // EXCA=1, voltage on cap = +/-VDD 
  ad7746.writeExcSetupRegister(0x2B); // EXCB=1, EXCA=1, voltage on cap = +/-VDD 
  //ad7746.writeConfigurationRegister(0xA1); // 62.1ms - VT conversion time; 62ms cap converstion time; continuous conversion
  ad7746.writeConfigurationRegister(0x01); // 20.1ms - VT conversion time; 11ms cap converstion time; continuous conversion
  //ad7746.writeCapDacARegister(0x00);
  ad7746.writeCapDacARegister(0x00); //MSB-enables DACA, 7b DAC; play around with CapDACA / or just calibrate ruler / or disconnect ruler and then calibrate
  //ad7746.writeCapDacARegister(0x37);

  // added calibration (+/-1pF; anything outside use CAPDAC)
  /*
  ad7746.write_register(0x0D,0x36); //LSB=32aF, 0.45pF is 14062 (0x36EE)
  ad7746.write_register(0x0E,0xEE); //LSB=32aF, 0.45pF is 14062 (0x36EE)
  */
  /*
  ad7746.write_register(0x0D,0xC9); //LSB=32aF, -0.45pF is -14062 (0xC912)
  ad7746.write_register(0x0E,0x12); //LSB=32aF, -0.45pF is -14062 (0xC912)
  */
  /*
  ad7746.write_register(0x0D,0xCF); //LSB=32aF, -0.4pF is -12500 (0xCF2C)
  ad7746.write_register(0x0E,0x2C); //LSB=32aF, -0.4pF is -12500 (0xCF2C)
  */
  Serial.println("AD7746 setup ok.");

}

void loop() {

  // choose CH1
  ad7746.writeCapSetupRegister(0x80); // CAPEN=1, CN2=0, DIFF=0, rest=0
  // Get capacitance measurement
  uint32_t capacitance_x = ad7746.getCapacitance();
  
  //Serial.println(capacitance);
  // Convert capacitance to farads (F) using calibration coefficient
  // capacitance code= 2^24-1; actual capacitance is 4, factor 8 is thus needed 
  double newValue_x = (((8.0/16777215.00) * capacitance_x) - 4); // capacitance is 0, capvalue is -4pF, 2^24 - 1 = 16777215
  //Serial.println(capacitance);
  // Print capacitance measurement in farads
  //Serial.print("Capacitance: ");
  Serial.println(newValue_x, 8); // Print with 8 decimal places for precision
  //Serial.println(" pF");

  // choose CH2
  ad7746.writeCapSetupRegister(0xC0); // CAPEN=1, CN2=1, DIFF=0, rest=0
  uint32_t capacitance_y = ad7746.getCapacitance();
  
  //Serial.println(capacitance);
  // Convert capacitance to farads (F) using calibration coefficient
  // capacitance code= 2^24-1; actual capacitance is 4, factor 8 is thus needed 
  double newValue_y = (((8.0/16777215.00) * capacitance_y) - 4); // capacitance is 0, capvalue is -4pF, 2^24 - 1 = 16777215
  //Serial.println(capacitance);
  // Print capacitance measurement in farads
  //Serial.print("Capacitance: ");
  Serial.println(newValue_y, 8); // Print with 8 decimal places for precision
  //Serial.println(" pF");

  //delay(100); // Delay for 1 second
}