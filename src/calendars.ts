/**
 *       JavaScript functions for the Fourmilab Calendar Converter
 *
 *                  by John Walker  --  September, MIM
 *              http://www.fourmilab.ch/documents/calendar/
 *
 *                This program is in the public domain.
 */

import { amod, deltat, equationOfTime, equinox, jwday, mod } from "./astro.js";
import {
  FRENCH_REVOLUTIONARY_EPOCH,
  GREGORIAN_EPOCH,
  HEBREW_EPOCH,
  ISLAMIC_EPOCH,
  MAYAN_COUNT_EPOCH,
  PERSIAN_EPOCH,
  TROPICAL_YEAR,
} from "./constants/index.js";

/**
 * Return Julian date of given weekday (0 = Sunday)
 * in the seven days ending on jd.
 *
 * @param weekday
 * @param jd
 */
export function weekday_before(weekday: number, jd: number): number {
  return jd - jwday(jd - weekday);
}

/**
 * Determine the Julian date for thi given parameters.
 *
 * @param weekday - Day of week desired, 0 = Sunday
 * @param jd - Julian date to begin search
 * @param direction - 1 = next weekday, -1 = last weekday
 * @param offset - Offset from jd to begin search
 */
export function search_weekday(weekday: number, jd: number, direction: number, offset: number): number {
  return weekday_before(weekday, jd + (direction * offset));
}

/**
 * Utility weekday functions, just wrappers for search_weekday.
 *
 * @param weekday
 * @param jd
 */
export function nearest_weekday(weekday: number, jd: number): number {
  return search_weekday(weekday, jd, 1, 3);
}

/**
 *
 * @param weekday
 * @param jd
 */
export function next_weekday(weekday: number, jd: number): number {
  return search_weekday(weekday, jd, 1, 7);
}

/**
 *
 * @param weekday
 * @param jd
 */
export function next_or_current_weekday(weekday: number, jd: number): number {
  return search_weekday(weekday, jd, 1, 6);
}

/**
 *
 * @param weekday
 * @param jd
 */
export function previous_weekday(weekday: number, jd: number): number {
  return search_weekday(weekday, jd, -1, 1);
}

/**
 *
 * @param weekday
 * @param jd
 */
export function previous_or_current_weekday(weekday: number, jd: number): number {
  return search_weekday(weekday, jd, 1, 0);
}

/**
 * Is a given year in the Gregorian calendar a leap year?
 *
 * @param year
 */
export function leap_gregorian(year: number): boolean {
  return ((year % 4) == 0) &&
    (!(((year % 100) == 0) && ((year % 400) != 0)));
}

/**
 * Determine Julian day number from Gregorian calendar date.
 *
 * @param year
 * @param month
 * @param day
 */
export function gregorian_to_jd(year: number, month: number, day: number): number {
  return (GREGORIAN_EPOCH - 1) +
    (365 * (year - 1)) +
    Math.floor((year - 1) / 4) +
    (-Math.floor((year - 1) / 100)) +
    Math.floor((year - 1) / 400) +
    Math.floor((((367 * month) - 362) / 12) +
      ((month <= 2) ? 0 :
        (leap_gregorian(year) ? -1 : -2)
      ) +
      day);
}

/**
 * Calculate Gregorian calendar date from Julian day.
 *
 * @param jd
 */
export function jd_to_gregorian(jd: number): [number, number, number] {
  const wjd = Math.floor(jd - 0.5) + 0.5;
  const depoch = wjd - GREGORIAN_EPOCH;
  const quadricent = Math.floor(depoch / 146097);
  const dqc = mod(depoch, 146097);
  const cent = Math.floor(dqc / 36524);
  const dcent = mod(dqc, 36524);
  const quad = Math.floor(dcent / 1461);
  const dquad = mod(dcent, 1461);
  const yindex = Math.floor(dquad / 365);

  let year = (quadricent * 400) + (cent * 100) + (quad * 4) + yindex;

  if (!((cent == 4) || (yindex == 4))) {
    year++;
  }
  const yearday = wjd - gregorian_to_jd(year, 1, 1);
  const leapadj = ((wjd < gregorian_to_jd(year, 3, 1)) ? 0
    :
    (leap_gregorian(year) ? 1 : 2)
  );
  const month = Math.floor((((yearday + leapadj) * 12) + 373) / 367);
  const day = (wjd - gregorian_to_jd(year, month, 1)) + 1;

  return [year, month, day];
}

/**
 * Return Julian day of given ISO year, week, and day.
 *
 * @param weekday
 * @param jd
 * @param nthweek
 */
export function n_weeks(weekday: number, jd: number, nthweek: number): number {
  let j = 7 * nthweek;

  if (nthweek > 0) {
    j += previous_weekday(weekday, jd);
  } else {
    j += next_weekday(weekday, jd);
  }
  return j;
}

/**
 *
 * @param year
 * @param week
 * @param day
 */
export function iso_to_julian(year: number, week: number, day: number): number {
  return day + n_weeks(0, gregorian_to_jd(year - 1, 12, 28), week);
}

/**
 * Return array of ISO (year, week, day) for Julian day.
 *
 * @param jd
 */
export function jd_to_iso(jd: number): [number, number, number] {
  let year = jd_to_gregorian(jd - 3)[0];
  if (jd >= iso_to_julian(year + 1, 1, 1)) {
    year++;
  }
  const week = Math.floor((jd - iso_to_julian(year, 1, 1)) / 7) + 1;
  let day = jwday(jd);
  if (day == 0) {
    day = 7;
  }
  return [year, week, day];
}

/**
 * Return Julian day of given ISO year, and day of year.
 *
 * @param year
 * @param day
 */
export function iso_day_to_julian(year: number, day: number): number {
  return (day - 1) + gregorian_to_jd(year, 1, 1);
}

/**
 * Return array of ISO (year, day_of_year) for Julian day.
 *
 * @param jd
 */
export function jd_to_iso_day(jd: number): [number, number] {
  const year = jd_to_gregorian(jd)[0];
  const day = Math.floor(jd - gregorian_to_jd(year, 1, 1)) + 1;

  return [year, day];
}

/**
 * Determine Julian day number from Julian calendar date
 *
 * @param year
 */
export function leap_julian(year: number): boolean {
  return mod(year, 4) == ((year > 0) ? 0 : 3);
}

/**
 *
 * @param year
 * @param month
 * @param day
 */
export function julian_to_jd(year: number, month: number, day: number): number {

  /* Adjust negative common era years to the zero-based notation we use.  */

  if (year < 1) {
    year++;
  }

  /* Algorithm as given in Meeus, Astronomical Algorithms, Chapter 7, page 61 */

  if (month <= 2) {
    year--;
    month += 12;
  }

  return ((Math.floor((365.25 * (year + 4716))) +
    Math.floor((30.6001 * (month + 1))) +
    day) - 1524.5);
}

/**
 * Calculate Julian calendar date from Julian day.
 *
 * @param td
 */
export function jd_to_julian(td: number): [number, number, number] {
  td += 0.5;

  const z = Math.floor(td);

  const a = z;
  const b = a + 1524;
  const c = Math.floor((b - 122.1) / 365.25);
  const d = Math.floor(365.25 * c);
  const e = Math.floor((b - d) / 30.6001);

  const month = Math.floor((e < 14) ? (e - 1) : (e - 13));
  let year = Math.floor((month > 2) ? (c - 4716) : (c - 4715));
  const day = b - d - Math.floor(30.6001 * e);

  /*  If year is less than 1, subtract one to convert from
      a zero based date system to the common era system in
      which the year -1 (1 B.C.E) is followed by year 1 (1 C.E.).  */

  if (year < 1) {
    year--;
  }

  return [year, month, day];
}

/**
 * Is a given Hebrew year a leap year?
 *
 * @param year
 */
export function hebrew_leap(year: number): boolean {
  return mod(((year * 7) + 1), 19) < 7;
}

/**
 * How many months are there in a Hebrew year (12 = normal, 13 = leap)
 *
 * @param year
 */
export function hebrew_year_months(year: number): number {
  return hebrew_leap(year) ? 13 : 12;
}

/**
 * Test for delay of start of new year and to avoid
 * Sunday, Wednesday, and Friday as start of the new year.
 *
 * @param year
 */
export function hebrew_delay_1(year: number): number {
  const months = Math.floor(((235 * year) - 234) / 19);
  const parts = 12084 + (13753 * months);
  let day = (months * 29) + Math.floor(parts / 25920);

  if (mod((3 * (day + 1)), 7) < 3) {
    day++;
  }
  return day;
}

/**
 * Check for delay in start of new year due to length of adjacent years
 *
 * @param year
 */
export function hebrew_delay_2(year: number): number {
  const last = hebrew_delay_1(year - 1);
  const present = hebrew_delay_1(year);
  const next = hebrew_delay_1(year + 1);

  return ((next - present) == 356) ? 2 :
    (((present - last) == 382) ? 1 : 0);
}

/**
 * How many days are in a Hebrew year?
 *
 * @param year
 */
export function hebrew_year_days(year: number): number {
  return hebrew_to_jd(year + 1, 7, 1) - hebrew_to_jd(year, 7, 1);
}

/**
 * How many days are in a given month of a given year
 *
 * @param year
 * @param month
 */
export function hebrew_month_days(year: number, month: number): number {
  //  First of all, dispose of fixed-length 29 day months

  if (month == 2 || month == 4 || month == 6 ||
    month == 10 || month == 13) {
    return 29;
  }

  //  If it's not a leap year, Adar has 29 days

  if (month == 12 && !hebrew_leap(year)) {
    return 29;
  }

  //  If it's Heshvan, days depend on length of year

  if (month == 8 && !(mod(hebrew_year_days(year), 10) == 5)) {
    return 29;
  }

  //  Similarly, Kislev varies with the length of year

  if (month == 9 && (mod(hebrew_year_days(year), 10) == 3)) {
    return 29;
  }

  //  Nope, it's a 30 day month

  return 30;
}

/**
 * Determine Julian day from Hebrew date.
 *
 * @param year
 * @param month
 * @param day
 */
export function hebrew_to_jd(year: number, month: number, day: number): number {
  const months = hebrew_year_months(year);
  let jd = HEBREW_EPOCH + hebrew_delay_1(year) +
    hebrew_delay_2(year) + day + 1;

  if (month < 7) {
    for (let mon = 7; mon <= months; mon++) {
      jd += hebrew_month_days(year, mon);
    }
    for (let mon = 1; mon < month; mon++) {
      jd += hebrew_month_days(year, mon);
    }
  } else {
    for (let mon = 7; mon < month; mon++) {
      jd += hebrew_month_days(year, mon);
    }
  }

  return jd;
}

/**
 * Convert Julian date to Hebrew date
 * This works by making multiple calls to
 * the inverse function, and is this very slow.
 *
 * @param jd
 */
export function jd_to_hebrew(jd: number): [number, number, number] {
  jd = Math.floor(jd) + 0.5;

  const count = Math.floor(((jd - HEBREW_EPOCH) * 98496.0) / 35975351.0);
  let year = count - 1;
  for (let i = count; jd >= hebrew_to_jd(i, 7, 1); i++) {
    year++;
  }
  const first = (jd < hebrew_to_jd(year, 1, 1)) ? 7 : 1;
  let month = first;
  for (let i = first; jd > hebrew_to_jd(year, i, hebrew_month_days(year, i)); i++) {
    month++;
  }
  const day = (jd - hebrew_to_jd(year, month, 1)) + 1;
  return [year, month, day];
}

/**
 * Determine Julian day and fraction of the
 * September equinox at the Paris meridian in
 * a given Gregorian year.
 *
 * @param year
 */
export function equinoxe_a_paris(year: number): number {
  //  September equinox in dynamical time
  const equJED = equinox(year, 2);

  //  Correct for delta T to obtain Universal time
  const equJD = equJED - (deltat(year) / (24 * 60 * 60));

  //  Apply the equation of time to yield the apparent time at Greenwich
  const equAPP = equJD + equationOfTime(equJED);

  /*  Finally, we must correct for the constant difference between
      the Greenwich meridian and that of Paris, 2�20'15" to the
      East.  */

  const dtParis = (2 + (20 / 60.0) + (15 / (60 * 60.0))) / 360;
  const equParis = equAPP + dtParis;

  return equParis;
}

/**
 * Calculate Julian day during which the
 * September equinox, reckoned from the Paris
 * meridian, occurred for a given Gregorian year.
 *
 * @param year
 */
export function paris_equinoxe_jd(year: number): number {
  const ep = equinoxe_a_paris(year);
  const epg = Math.floor(ep - 0.5) + 0.5;

  return epg;
}

/**
 * Determine the year in the French
 * revolutionary calendar in which a
 * given Julian day falls.  Returns an
 * array of two elements:
 * [0]  Ann�e de la R�volution
 * [1]  Julian day number containing equinox for this year.
 *
 * @param jd
 */
export function annee_da_la_revolution(jd: number): [number, number] {
  let guess = jd_to_gregorian(jd)[0] - 2,
    lasteq, nexteq;

  lasteq = paris_equinoxe_jd(guess);
  while (lasteq > jd) {
    guess--;
    lasteq = paris_equinoxe_jd(guess);
  }
  nexteq = lasteq - 1;
  while (!((lasteq <= jd) && (jd < nexteq))) {
    lasteq = nexteq;
    guess++;
    nexteq = paris_equinoxe_jd(guess);
  }
  const adr = Math.round((lasteq - FRENCH_REVOLUTIONARY_EPOCH) / TROPICAL_YEAR) + 1;

  return [adr, lasteq];
}

/**
 * Calculate date in the French Revolutionary
 * calendar from Julian day.  The five or six
 * "sansculottides" are considered a thirteenth
 * month in the results of this function.
 *
 * @param jd
 */
export function jd_to_french_revolutionary(jd: number): [number, number, number, number] {
  jd = Math.floor(jd) + 0.5;

  const adr = annee_da_la_revolution(jd);
  const an = adr[0];
  const equinoxe = adr[1];
  const mois = Math.floor((jd - equinoxe) / 30) + 1;
  let jour = (jd - equinoxe) % 30;
  const decade = Math.floor(jour / 10) + 1;
  jour = (jour % 10) + 1;

  return [an, mois, decade, jour];
}

/**
 * Obtain Julian day from a given French
 * Revolutionary calendar date.
 *
 * @param an
 * @param mois
 * @param decade
 * @param jour
 */
export function french_revolutionary_to_jd(an: number, mois: number, decade: number, jour: number): number {
  let adr, guess;

  guess = FRENCH_REVOLUTIONARY_EPOCH + (TROPICAL_YEAR * ((an - 1) - 1));
  adr = [an - 1, 0];

  while (adr[0] < an) {
    adr = annee_da_la_revolution(guess);
    guess = adr[1] + (TROPICAL_YEAR + 2);
  }
  const equinoxe = adr[1];

  const jd = equinoxe + (30 * (mois - 1)) + (10 * (decade - 1)) + (jour - 1);
  return jd;
}

/**
 * Is a given year a leap year in the Islamic calendar?
 *
 * @param year
 */
export function leap_islamic(year: number): boolean {
  return (((year * 11) + 14) % 30) < 11;
}

/**
 * Determine Julian day from Islamic date
 *
 * @param year
 * @param month
 * @param day
 */
export function islamic_to_jd(year: number, month: number, day: number): number {
  return (day +
    Math.ceil(29.5 * (month - 1)) +
    (year - 1) * 354 +
    Math.floor((3 + (11 * year)) / 30) +
    ISLAMIC_EPOCH) - 1;
}

/**
 * Calculate Islamic date from Julian day
 *
 * @param jd
 */
export function jd_to_islamic(jd: number): [number, number, number] {
  jd = Math.floor(jd) + 0.5;

  const year = Math.floor(((30 * (jd - ISLAMIC_EPOCH)) + 10646) / 10631);
  const month = Math.min(12,
    Math.ceil((jd - (29 + islamic_to_jd(year, 1, 1))) / 29.5) + 1);
  const day = (jd - islamic_to_jd(year, month, 1)) + 1;

  return [year, month, day];
}

/**
 * Determine Julian day and fraction of the
 * March equinox at the Tehran meridian in
 * a given Gregorian year.
 *
 * @param year
 */
export function tehran_equinox(year: number): number {
  //  March equinox in dynamical time
  const equJED = equinox(year, 0);

  //  Correct for delta T to obtain Universal time
  const equJD = equJED - (deltat(year) / (24 * 60 * 60));

  //  Apply the equation of time to yield the apparent time at Greenwich
  const equAPP = equJD + equationOfTime(equJED);

  /*  Finally, we must correct for the constant difference between
      the Greenwich meridian andthe time zone standard for
Iran Standard time, 52�30' to the East.  */

  const dtTehran = (52 + (30 / 60.0) + (0 / (60.0 * 60.0))) / 360;
  const equTehran = equAPP + dtTehran;

  return equTehran;
}

/**
 * Calculate Julian day during which the
 * March equinox, reckoned from the Tehran
 * meridian, occurred for a given Gregorian year.
 *
 * @param year
 */
export function tehran_equinox_jd(year: number): number {
  const ep = tehran_equinox(year);
  const epg = Math.floor(ep);

  return epg;
}

/**
 * Determine the year in the Persian
 * astronomical calendar in which a
 * given Julian day falls.  Returns an
 * array of two elements:
 * [0]  Persian year
 * [1]  Julian day number containing equinox for this year.
 *
 * @param jd
 */
export function persiana_year(jd: number): [number, number] {
  let guess = jd_to_gregorian(jd)[0] - 2,
    lasteq, nexteq;

  lasteq = tehran_equinox_jd(guess);
  while (lasteq > jd) {
    guess--;
    lasteq = tehran_equinox_jd(guess);
  }
  nexteq = lasteq - 1;
  while (!((lasteq <= jd) && (jd < nexteq))) {
    lasteq = nexteq;
    guess++;
    nexteq = tehran_equinox_jd(guess);
  }
  const adr = Math.round((lasteq - PERSIAN_EPOCH) / TROPICAL_YEAR) + 1;

  return [adr, lasteq];
}

/**
 * Calculate date in the Persian astronomical
 * calendar from Julian day.
 *
 * @param jd
 */
export function jd_to_persiana(jd: number): [number, number, number] {
  jd = Math.floor(jd) + 0.5;

  const adr = persiana_year(jd);
  const year = adr[0];
  // const equinox = adr[1];
  // day = Math.floor((jd - equinox) / 30) + 1;

  const yday = (Math.floor(jd) - persiana_to_jd(year, 1, 1)) + 1;
  const month = (yday <= 186) ? Math.ceil(yday / 31) : Math.ceil((yday - 6) / 30);
  const day = (Math.floor(jd) - persiana_to_jd(year, month, 1)) + 1;

  return [year, month, day];
}

/**
 * Obtain Julian day from a given Persian astronomical calendar date.
 *
 * @param year
 * @param month
 * @param day
 */
export function persiana_to_jd(year: number, month: number, day: number): number {
  let adr, guess;

  guess = (PERSIAN_EPOCH - 1) + (TROPICAL_YEAR * ((year - 1) - 1));
  adr = [year - 1, 0];

  while (adr[0] < year) {
    adr = persiana_year(guess);
    guess = adr[1] + (TROPICAL_YEAR + 2);
  }
  const equinox = adr[1];

  const jd = equinox +
    ((month <= 7) ?
      ((month - 1) * 31) :
      (((month - 1) * 30) + 6)
    ) +
    (day - 1);
  return jd;
}

/**
 * Is a given year a leap year in the Persian astronomical calendar?
 *
 * @param year
 */
export function leap_persiana(year: number): boolean {
  return (persiana_to_jd(year + 1, 1, 1) -
    persiana_to_jd(year, 1, 1)) > 365;
}

/**
 * Is a given year a leap year in the Persian calendar?
 *
 * @param year
 */
export function leap_persian(year: number): boolean {
  return ((((((year - ((year > 0) ? 474 : 473)) % 2820) + 474) + 38) * 682) % 2816) < 682;
}

/**
 * Determine Julian day from Persian date
 *
 * @param year
 * @param month
 * @param day
 */
export function persian_to_jd(year: number, month: number, day: number): number {
  const epbase = year - ((year >= 0) ? 474 : 473);
  const epyear = 474 + mod(epbase, 2820);

  return day +
    ((month <= 7) ?
      ((month - 1) * 31) :
      (((month - 1) * 30) + 6)
    ) +
    Math.floor(((epyear * 682) - 110) / 2816) +
    (epyear - 1) * 365 +
    Math.floor(epbase / 2820) * 1029983 +
    (PERSIAN_EPOCH - 1);
}

/**
 * Calculate Persian date from Julian day
 *
 * @param jd
 */
export function jd_to_persian(jd: number): [number, number, number] {
  let ycycle, aux1, aux2;

  jd = Math.floor(jd) + 0.5;

  const depoch = jd - persian_to_jd(475, 1, 1);
  const cycle = Math.floor(depoch / 1029983);
  const cyear = mod(depoch, 1029983);
  if (cyear == 1029982) {
    ycycle = 2820;
  } else {
    aux1 = Math.floor(cyear / 366);
    aux2 = mod(cyear, 366);
    ycycle = Math.floor(((2134 * aux1) + (2816 * aux2) + 2815) / 1028522) +
      aux1 + 1;
  }
  let year = ycycle + (2820 * cycle) + 474;
  if (year <= 0) {
    year--;
  }
  const yday = (jd - persian_to_jd(year, 1, 1)) + 1;
  const month = (yday <= 186) ? Math.ceil(yday / 31) : Math.ceil((yday - 6) / 30);
  const day = (jd - persian_to_jd(year, month, 1)) + 1;
  return [year, month, day];
}

/**
 * Determine Julian day from Mayan long count
 *
 * @param baktun
 * @param katun
 * @param tun
 * @param uinal
 * @param kin
 */
export function mayan_count_to_jd(baktun: number, katun: number, tun: number, uinal: number, kin: number): number {
  return MAYAN_COUNT_EPOCH +
    (baktun * 144000) +
    (katun * 7200) +
    (tun * 360) +
    (uinal * 20) +
    kin;
}

/**
 * Calculate Mayan long count from Julian day
 *
 * @param jd
 */
export function jd_to_mayan_count(jd: number): [number, number, number, number, number] {
  let d;

  jd = Math.floor(jd) + 0.5;
  d = jd - MAYAN_COUNT_EPOCH;
  const baktun = Math.floor(d / 144000);
  d = mod(d, 144000);
  const katun = Math.floor(d / 7200);
  d = mod(d, 7200);
  const tun = Math.floor(d / 360);
  d = mod(d, 360);
  const uinal = Math.floor(d / 20);
  const kin = mod(d, 20);

  return [baktun, katun, tun, uinal, kin];
}

/**
 * Determine Mayan Haab "month" and day from Julian day
 *
 * @param jd
 */
export function jd_to_mayan_haab(jd: number): [number, number] {
  jd = Math.floor(jd) + 0.5;

  const lcount = jd - MAYAN_COUNT_EPOCH;
  const day = mod(lcount + 8 + ((18 - 1) * 20), 365);

  return [Math.floor(day / 20) + 1, mod(day, 20)];
}

/**
 * Determine Mayan Tzolkin "month" and day from Julian day
 *
 * @param jd
 */
export function jd_to_mayan_tzolkin(jd: number): [number, number] {
  jd = Math.floor(jd) + 0.5;

  const lcount = jd - MAYAN_COUNT_EPOCH;

  return [amod(lcount + 20, 20), amod(lcount + 4, 13)];
}

/**
 * Obtain Julian day for Indian Civil date
 *
 * @param year
 * @param month
 * @param day
 */
export function indian_civil_to_jd(year: number, month: number, day: number): number {
  let jd, m;

  const gyear = year + 78;
  const leap = leap_gregorian(gyear);     // Is this a leap year ?
  const start = gregorian_to_jd(gyear, 3, leap ? 21 : 22);
  const Caitra = leap ? 31 : 30;

  if (month == 1) {
    jd = start + (day - 1);
  } else {
    jd = start + Caitra;
    m = month - 2;
    m = Math.min(m, 5);
    jd += m * 31;
    if (month >= 8) {
      m = month - 7;
      jd += m * 30;
    }
    jd += day - 1;
  }

  return jd;
}

/**
 * Calculate Indian Civil date from Julian day
 *
 * @param jd
 */
export function jd_to_indian_civil(jd: number): [number, number, number] {
  let year, yday, mday, month, day;

  const Saka = 79 - 1;                    // Offset in years from Saka era to Gregorian epoch
  const start = 80;                       // Day offset between Saka and Gregorian

  jd = Math.floor(jd) + 0.5;
  const greg = jd_to_gregorian(jd);       // Gregorian date for Julian day
  const leap = leap_gregorian(greg[0]);   // Is this a leap year?
  year = greg[0] - Saka;            // Tentative year in Saka era
  const greg0 = gregorian_to_jd(greg[0], 1, 1); // JD at start of Gregorian year
  yday = jd - greg0;                // Day number (0 based) in Gregorian year
  const Caitra = leap ? 31 : 30;          // Days in Caitra this year

  if (yday < start) {
    //  Day is at the end of the preceding Saka year
    year--;
    yday += Caitra + (31 * 5) + (30 * 3) + 10 + start;
  }

  yday -= start;
  if (yday < Caitra) {
    month = 1;
    day = yday + 1;
  } else {
    mday = yday - Caitra;
    if (mday < (31 * 5)) {
      month = Math.floor(mday / 31) + 2;
      day = (mday % 31) + 1;
    } else {
      mday -= 31 * 5;
      month = Math.floor(mday / 30) + 7;
      day = (mday % 30) + 1;
    }
  }

  return [year, month, day];
}
