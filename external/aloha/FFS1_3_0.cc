// This File is Automatically generated by ALOHA
// 

#include "FFS1_3_0.h"
#include "FFS1_0.h"
#include "FFS3_0.h"

void FFS1_3_0(std::complex<double> F1[], std::complex<double> F2[], std::complex<double> S3[],
    std::complex<double> COUP1, std::complex<double> COUP2, std::complex<double> & vertex)
{
  std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> tmp; 
  FFS1_0(F1, F2, S3, COUP1, vertex); 
  FFS3_0(F1, F2, S3, COUP2, tmp); 
  vertex = vertex + tmp; 
}

