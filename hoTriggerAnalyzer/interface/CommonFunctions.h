#ifndef __HOMUON__COMMONFUNCTIONS_H__
#define __HOMUON__COMMONFUNCTIONS_H__

/*
 * Common Functions Class
 * Author Chris Anelli
 * 6.13.2013
 */


class CommonFunctions {  
  
 public:

  /*
   * Takes the difference of two phis, makes
   * sure they are not more than 2 pi.
   */
  float WrapCheck(float phi1, float phi2);

 private:

  
};

#endif
