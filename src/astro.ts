/**
 *            JavaScript functions for positional astronomy
 *
 *                  by John Walker  --  September, MIM
 *                       http://www.fourmilab.ch/
 *
 *                This program is in the public domain.
 */

import {
  DEGREES_180,
  DEGREES_360,
  deltaTtab,
  EquinoxpTerms,
  J2000,
  JDE0tab1000,
  JDE0tab2000,
  JULIAN_CENTURY,
  JULIAN_MILLENNIUM,
  nutArgCoeff,
  nutArgMult,
  oterms,
  PI,
  PI2,
  SECONDS_PER_HOUR,
} from "./constants";

/**
 * Arc-seconds to radians.
 *
 * @param seconds
 */
export function astor(seconds: number): number {
  return seconds * (PI / (DEGREES_180 * SECONDS_PER_HOUR));
}

/**
 * Degrees to radians.
 *
 * @param degrees
 */
export function dtr(degrees: number): number {
  return (degrees * PI) / DEGREES_180;
}

/**
 * Radians to degrees.
 *
 * @param radians
 */
export function rtd(radians: number): number {
  return (radians * DEGREES_180) / PI;
}

/**
 * Range reduce angle in degrees.
 *
 * @param angle
 */
export function fixangle(angle: number): number {
  return angle - DEGREES_360 * (Math.floor(angle / DEGREES_360));
}

/**
 * Range reduce angle in radians.
 *
 * @param angle
 */
export function fixangr(angle: number): number {
  return angle - PI2 * (Math.floor(angle / PI2));
}

/**
 * Sine of an angle in degrees.
 *
 * @param degrees
 */
export function dsin(degrees: number): number {
  return Math.sin(dtr(degrees));
}

/**
 * Cosine of an angle in degrees.
 *
 * @param degrees
 */
export function dcos(degrees: number): number {
  return Math.cos(dtr(degrees));
}

/**
 * Modulus function which works for non-integers.
 *
 * @param a
 * @param b
 */
export function mod(a: number, b: number): number {
  return a - (b * Math.floor(a / b));
}

/**
 * Modulus function which returns numerator if modulus is zero.
 *
 * @param a
 * @param b
 */
export function amod(a: number, b: number): number {
  return mod(a - 1, b) + 1;
}

/**
 * Convert Julian time to hour, minutes, and seconds,
 * returned as a three-element array.
 *
 * @param julian
 */
export function jhms(julian: number): [number, number, number] {
  /* Astronomical to civil */
  julian += 0.5;

  const ij = ((julian - Math.floor(julian)) * 86400.0) + 0.5;

  return [
    Math.floor(ij / 3600),
    Math.floor((ij / 60) % 60),
    Math.floor(ij % 60),
  ];
}

/**
 * Calculate day of week from Julian day.
 *
 * @param julian
 */
export function jwday(julian: number): number {
  return mod(Math.floor((julian + 1.5)), 7);
}

/**
 * Calculate the obliquity of the ecliptic for a given
 * Julian date.  This uses Laskar's tenth-degree
 * polynomial fit (J. Laskar, Astronomy and
 * Astrophysics, Vol. 157, page 68 [1986]) which is
 * accurate to within 0.01 arc second between AD 1000
 * and AD 3000, and within a few seconds of arc for
 * +/-10000 years around AD 2000.  If we're outside the
 * range in which this fit is valid (deep time) we
 * simply return the J2000 value of the obliquity, which
 * happens to be almost precisely the mean.
 *
 * @param julian
 */
export function obliqeq(julian: number): number {
  const u = (julian - J2000) / (JULIAN_CENTURY * 100);
  let v = u;
  let eps = 23 + (26 / 60.0) + (21.448 / 3600.0);

  if (Math.abs(u) >= 1.0) return eps;

  for (let i = 0; i < 10; i++) {
    eps += (oterms[i] / 3600.0) * v;

    v *= u;
  }

  return eps;
}

/**
 * Calculate the nutation in longitude, deltaPsi, and
 * obliquity, deltaEpsilon for a given Julian date
 * jd.  Results are returned as a two element Array
 * giving (deltaPsi, deltaEpsilon) in degrees.
 *
 * @param julian
 */
export function nutation(julian: number): [number, number] {
  let dp = 0, de = 0, ang;

  const t = (julian - 2451545.0) / 36525.0;
  const ta = [];
  const t2 = t * t;
  const t3 = t * t2;

  /* Calculate angles.  The correspondence between the elements
     of our array and the terms cited in Meeus are:

     ta[0] = D  ta[0] = M  ta[2] = M'  ta[3] = F  ta[4] = \Omega

  */

  ta[0] = dtr(297.850363 + 445267.11148 * t - 0.0019142 * t2 +
    t3 / 189474.0);
  ta[1] = dtr(357.52772 + 35999.05034 * t - 0.0001603 * t2 -
    t3 / 300000.0);
  ta[2] = dtr(134.96298 + 477198.867398 * t + 0.0086972 * t2 +
    t3 / 56250.0);
  ta[3] = dtr(93.27191 + 483202.017538 * t - 0.0036825 * t2 +
    t3 / 327270);
  ta[4] = dtr(125.04452 - 1934.136261 * t + 0.0020708 * t2 +
    t3 / 450000.0);

  /* Range reduce the angles in case the sine and cosine functions
     don't do it as accurately or quickly. */

  for (let i = 0; i < 5; i++) {
    ta[i] = fixangr(ta[i]);
  }

  const to10 = t / 10.0;

  for (let i = 0; i < 63; i++) {
    ang = 0;

    for (let j = 0; j < 5; j++) {
      if (nutArgMult[(i * 5) + j] != 0) {
        ang += nutArgMult[(i * 5) + j] * ta[j];
      }
    }

    dp += (nutArgCoeff[(i * 4) + 0] + nutArgCoeff[(i * 4) + 1] * to10) * Math.sin(ang);

    de += (nutArgCoeff[(i * 4) + 2] + nutArgCoeff[(i * 4) + 3] * to10) * Math.cos(ang);
  }

  /* Return the result, converting from ten thousandths of arc
     seconds to radians in the process. */

  const deltaPsi = dp / (3600.0 * 10000.0);
  const deltaEpsilon = de / (3600.0 * 10000.0);

  return [deltaPsi, deltaEpsilon];
}

/**
 * Convert celestial (ecliptical) longitude and
 * latitude into right ascension (in degrees) and
 * declination.  We must supply the time of the
 * conversion in order to compensate correctly for the
 * varying obliquity of the ecliptic over time.
 * The right ascension and declination are returned
 * as a two-element Array in that order.
 *
 * @param jd
 * @param Lambda
 * @param Beta
 */
export function ecliptoeq(jd: number, Lambda: number, Beta: number): [number, number] {
  /* Obliquity of the ecliptic. */

  const eps = dtr(obliqeq(jd));

  // Ra = rtd(Math.atan2((Math.cos(eps) * Math.sin(dtr(Lambda)) -
  //     (Math.tan(dtr(Beta)) * Math.sin(eps))),
  // Math.cos(dtr(Lambda))));
  const Ra = fixangle(rtd(Math.atan2((Math.cos(eps) * Math.sin(dtr(Lambda)) -
      (Math.tan(dtr(Beta)) * Math.sin(eps))),
  Math.cos(dtr(Lambda)))));
  const Dec = rtd(Math.asin((Math.sin(eps) * Math.sin(dtr(Lambda)) * Math.cos(dtr(Beta))) +
    (Math.sin(dtr(Beta)) * Math.cos(eps))));

  return [Ra, Dec];
}

/**
 * Determine the difference, in seconds, between
 * Dynamical time and Universal time.
 *
 * @param year
 */
export function deltat(year: number): number {
  let dt, f, i, t;

  if ((year >= 1620) && (year <= 2000)) {
    i = Math.floor((year - 1620) / 2);

    f = ((year - 1620) / 2) - i;  /* Fractional part of year */

    dt = deltaTtab[i] + ((deltaTtab[i + 1] - deltaTtab[i]) * f);
  } else {
    t = (year - 2000) / 100;

    if (year < 948) {
      dt = 2177 + (497 * t) + (44.1 * t * t);
    } else {
      dt = 102 + (102 * t) + (25.3 * t * t);

      if ((year > 2000) && (year < 2100)) {
        dt += 0.37 * (year - 2100);
      }
    }
  }

  return dt;
}

/**
 * Determine the Julian Ephemeris Day of an equinox or solstice.
 * The "which" argument selects the item to be computed:
 * 0   March equinox
 * 1   June solstice
 * 2   September equinox
 * 3   December solstice
 *
 * @param year
 * @param which
 */
export function equinox(year: number, which: number): number {
  let i, j, JDE0tab, S, Y;

  /*  Initialise terms for mean equinox and solstices.  We
      have two sets: one for years prior to 1000 and a second
      for subsequent years.  */

  if (year < 1000) {
    JDE0tab = JDE0tab1000;
    Y = year / 1000;
  } else {
    JDE0tab = JDE0tab2000;
    Y = (year - 2000) / 1000;
  }

  const JDE0 = JDE0tab[which][0] +
    (JDE0tab[which][1] * Y) +
    (JDE0tab[which][2] * Y * Y) +
    (JDE0tab[which][3] * Y * Y * Y) +
    (JDE0tab[which][4] * Y * Y * Y * Y);

  //document.debug.log.value += "JDE0 = " + JDE0 + "\n";

  const T = (JDE0 - 2451545.0) / 36525;
  //document.debug.log.value += "T = " + T + "\n";
  const W = (35999.373 * T) - 2.47;
  //document.debug.log.value += "W = " + W + "\n";
  const deltaL = 1 + (0.0334 * dcos(W)) + (0.0007 * dcos(2 * W));
  //document.debug.log.value += "deltaL = " + deltaL + "\n";

  //  Sum the periodic terms for time T

  S = 0;
  for (i = j = 0; i < 24; i++) {
    S += EquinoxpTerms[j] * dcos(EquinoxpTerms[j + 1] + (EquinoxpTerms[j + 2] * T));
    j += 3;
  }

  //document.debug.log.value += "S = " + S + "\n";
  //document.debug.log.value += "Corr = " + ((S * 0.00001) / deltaL) + "\n";

  const JDE = JDE0 + ((S * 0.00001) / deltaL);

  return JDE;
}

/**
 * Position of the Sun.  Please see the comments
 * on the return statement at the end of this function
 * which describe the array it returns.  We return
 * intermediate values because they are useful in a
 * variety of other contexts.
 *
 * @param julian
 */
export function sunpos(julian: number): [
  number, number, number, number, number, number, number, number, number, number, number, number,
] {
  let L0, M, Alpha, AlphaApp;

  const T = (julian - J2000) / JULIAN_CENTURY;
  //document.debug.log.value += "Sunpos.  T = " + T + "\n";
  const T2 = T * T;
  L0 = 280.46646 + (36000.76983 * T) + (0.0003032 * T2);
  //document.debug.log.value += "L0 = " + L0 + "\n";
  L0 = fixangle(L0);
  //document.debug.log.value += "L0 = " + L0 + "\n";
  M = 357.52911 + (35999.05029 * T) + (-0.0001537 * T2);
  //document.debug.log.value += "M = " + M + "\n";
  M = fixangle(M);
  //document.debug.log.value += "M = " + M + "\n";
  const e = 0.016708634 + (-0.000042037 * T) + (-0.0000001267 * T2);
  //document.debug.log.value += "e = " + e + "\n";
  const C = ((1.914602 + (-0.004817 * T) + (-0.000014 * T2)) * dsin(M)) +
    ((0.019993 - (0.000101 * T)) * dsin(2 * M)) +
    (0.000289 * dsin(3 * M));
  //document.debug.log.value += "C = " + C + "\n";
  const sunLong = L0 + C;
  //document.debug.log.value += "sunLong = " + sunLong + "\n";
  const sunAnomaly = M + C;
  //document.debug.log.value += "sunAnomaly = " + sunAnomaly + "\n";
  const sunR = (1.000001018 * (1 - (e * e))) / (1 + (e * dcos(sunAnomaly)));
  //document.debug.log.value += "sunR = " + sunR + "\n";
  const Omega = 125.04 - (1934.136 * T);
  //document.debug.log.value += "Omega = " + Omega + "\n";
  const Lambda = sunLong + (-0.00569) + (-0.00478 * dsin(Omega));
  //document.debug.log.value += "Lambda = " + Lambda + "\n";
  const epsilon0 = obliqeq(julian);
  //document.debug.log.value += "epsilon0 = " + epsilon0 + "\n";
  const epsilon = epsilon0 + (0.00256 * dcos(Omega));
  //document.debug.log.value += "epsilon = " + epsilon + "\n";
  Alpha = rtd(Math.atan2(dcos(epsilon0) * dsin(sunLong), dcos(sunLong)));
  //document.debug.log.value += "Alpha = " + Alpha + "\n";
  Alpha = fixangle(Alpha);
  ////document.debug.log.value += "Alpha = " + Alpha + "\n";
  const Delta = rtd(Math.asin(dsin(epsilon0) * dsin(sunLong)));
  ////document.debug.log.value += "Delta = " + Delta + "\n";
  AlphaApp = rtd(Math.atan2(dcos(epsilon) * dsin(Lambda), dcos(Lambda)));
  //document.debug.log.value += "AlphaApp = " + AlphaApp + "\n";
  AlphaApp = fixangle(AlphaApp);
  //document.debug.log.value += "AlphaApp = " + AlphaApp + "\n";
  const DeltaApp = rtd(Math.asin(dsin(epsilon) * dsin(Lambda)));
  //document.debug.log.value += "DeltaApp = " + DeltaApp + "\n";

  return [                 //  Angular quantities are expressed in decimal degrees
    L0,                           //  [0] Geometric mean longitude of the Sun
    M,                            //  [1] Mean anomaly of the Sun
    e,                            //  [2] Eccentricity of the Earth's orbit
    C,                            //  [3] Sun's equation of the Centre
    sunLong,                      //  [4] Sun's true longitude
    sunAnomaly,                   //  [5] Sun's true anomaly
    sunR,                         //  [6] Sun's radius vector in AU
    Lambda,                       //  [7] Sun's apparent longitude at true equinox of the date
    Alpha,                        //  [8] Sun's true right ascension
    Delta,                        //  [9] Sun's true declination
    AlphaApp,                     // [10] Sun's apparent right ascension
    DeltaApp,                      // [11] Sun's apparent declination
  ];
}

/**
 * Compute equation of time for a given moment.
 * Returns the equation of time as a fraction of a day.
 *
 * @param julian
 */
export function equationOfTime(julian: number): number {
  const tau = (julian - J2000) / JULIAN_MILLENNIUM;

  let L0 = 280.4664567 + (360007.6982779 * tau) +
    (0.03032028 * tau * tau) +
    ((tau * tau * tau) / 49931) +
    (-((tau * tau * tau * tau) / 15300)) +
    (-((tau * tau * tau * tau * tau) / 2000000));

  L0 = fixangle(L0);

  const alpha = sunpos(julian)[10];

  const deltaPsi = nutation(julian)[0];

  const epsilon = obliqeq(julian) + nutation(julian)[1];

  let E = L0 + (-0.0057183) + (-alpha) + (deltaPsi * dcos(epsilon));

  E = E - 20.0 * (Math.floor(E / 20.0));

  E = E / (24 * 60);

  return E;
}
