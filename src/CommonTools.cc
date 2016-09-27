#include "CommonTools.h"

bool CommonTools::hardwareSimulation=true;

/* Function which simulate the HardWare representation of the values : manage UNSIGNED and SIGNED (2's complement) overflows and accuracy according to the available dynamic of the binary word */
double CommonTools::binning(double fNumber, int nMSBpowOfTwo, int nBits, HW_SIGN_TYPE signType)
{
  if (!hardwareSimulation)
    //If the Hardware binning simulation is not asked, return directly the original number
    return fNumber;

  if (signType == UNSIGNED && fNumber < 0)
    {
      //Bad interpretation, a negative number is stored in an UNSIGNED format (sign lost)
      fNumber = -fNumber;
    }
  
  int nLSBpowOfTwo;
	
  //Process the power of two of the LSB for the binary representation
  if (signType == UNSIGNED)
    {
      //If UNSIGNED
      nLSBpowOfTwo = nMSBpowOfTwo - (nBits-1);
    }	
  else
    {
      //If SIGNED, 1 bit is used for the sign
      nLSBpowOfTwo = nMSBpowOfTwo - (nBits-2);
    }

  /* Accuracy Simulation */

  //Divide the number by the power of two of the LSB => the integer part of the new number is the value we are looking for
  fNumber = fNumber / pow(2, nLSBpowOfTwo);
	
  //Remove the fractionnal part by rounding down (for both positive and negative values), this simulate the HW truncature
  fNumber = floor(fNumber);
	
  //Multiply the number by the power of two of the LSB to get the correct float value
  fNumber = fNumber * pow(2, nLSBpowOfTwo);


  double fBinnedNumber = fNumber;

  /* Overflow Simulation */

  if (signType == UNSIGNED)
    {
      //If the number is in UNSIGNED representation
      fNumber = fmod(fNumber, pow(2, nMSBpowOfTwo+1));
    }
  else
    {
      //If the number is in SIGNED representation (2's complement)
      
      double fTempResult = fNumber - pow(2, nMSBpowOfTwo+1); //substract the possible range to the number

      if (fTempResult >= 0)
        {
          //If there is an overflow, it's a positive one
          fNumber = fmod(fTempResult, pow(2, nMSBpowOfTwo+2)) - pow(2, nMSBpowOfTwo+1);
        }
      else
        {
          //If there is an overflow, it's a negative one (2's complement format has an asymetric range for positive and negative values)
          fNumber = fmod(fTempResult + pow(2, nLSBpowOfTwo), pow(2, nMSBpowOfTwo+2)) - pow(2, nLSBpowOfTwo) + pow(2, nMSBpowOfTwo+1);
        }
    }

  //If the new number is different from the previous one, an HW overflow occured
  if (fNumber != fBinnedNumber)
    {
      cout<<"WARNING HW overflow for the value : "<<fBinnedNumber<<" resulting value : "<<fNumber<<" (diff= "<<fBinnedNumber-fNumber<<")"<<endl;
    }
	
  return fNumber;
}
