#include "HoMuonTrigger/hoTriggerAnalyzer/interface/CommonFunctions.h"

/*
 * The Common Function Class contains  
 * functions that are useful in 
 * my different circumstances, and
 * that I have used often in my analyses. 
 *
 * Created by Christopher Anelli
 * On 6.15.2014
*/

#include "math.h"

/*
 * Wrap check calcuates the difference between two phi's,
 * making sure they are not more than 2 pi apart.
 */

float CommonFunctions::WrapCheck(float phi1, float phi2){
  //double M_PI = (double) 3.14;
  float delta_phi = phi1 - phi2;
  if(delta_phi < -M_PI){
    return (2*M_PI + delta_phi);
  }
  if(delta_phi > M_PI){
    return (delta_phi - 2*M_PI);
  }
  return delta_phi;
};
