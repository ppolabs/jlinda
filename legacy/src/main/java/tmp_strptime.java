public class GlobalMembersTmp_strptime
{


	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if ! HAVE_LOCALTIME_R && ! localtime_r
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if _LIBC
	//C++ TO JAVA CONVERTER TODO TASK: Alternate #define macros with the same name cannot be converted to Java:
	//#define localtime_r __localtime_r
	//#else
	// Approximate localtime_r as best we can in its absence.  
	//C++ TO JAVA CONVERTER TODO TASK: Alternate #define macros with the same name cannot be converted to Java:
	//#define localtime_r my_localtime_r
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro __P was defined in alternate ways and cannot be replaced in-line:
	static struct tm * my_localtime_r __P((const time_t *t, struct tm *tp))
	{
	  tm l = localtime (t);
	  if (l == null)
		return 0;
	  *tp = l;
	  return tp;
	}
	//#endif // ! _LIBC
	//#endif // ! HAVE_LOCALTIME_R && ! defined (localtime_r)

	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if WIN32
	  // Jia defined this for windows.
	  // Bert Kampes, 24-Aug-2005
	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//#define strncasecmp _strnicmp
	//#endif



	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//#define match_char(ch1, ch2) if (ch1 != ch2) return NULL
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if __GNUC__ && __GNUC__ >= 2
	//C++ TO JAVA CONVERTER TODO TASK: Alternate #define macros with the same name cannot be converted to Java:
	//#define match_string(cs1, s2) ({ size_t len = strlen (cs1); int result = strncasecmp ((cs1), (s2), len) == 0; if (result) (s2) += len; result; })
	//#else
	// Oh come on.  Get a reasonable compiler.  
	//C++ TO JAVA CONVERTER TODO TASK: Alternate #define macros with the same name cannot be converted to Java:
	//#define match_string(cs1, s2) (strncasecmp ((cs1), (s2), strlen (cs1)) ? 0 : ((s2) += strlen (cs1), 1))
	//#endif
	// We intentionally do not use isdigit() for testing because this will
	//   lead to problems with the wide character version.  
	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//#define get_number(from, to) do { val = 0; while (*rp == ' ') ++rp; if (*rp < '0' || *rp > '9') return NULL; do { val *= 10; val += *rp++ - '0'; } while (val * 10 <= to && *rp >= '0' && *rp <= '9'); if (val < from || val > to) return NULL; } while (0)
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if _NL_CURRENT
	//C++ TO JAVA CONVERTER TODO TASK: Alternate #define macros with the same name cannot be converted to Java:
	//#define get_alt_number (from, to) do { if (*decided != raw) { const char *alts = _NL_CURRENT (LC_TIME, ALT_DIGITS); val = 0; while (*alts != '\0') { size_t len = strlen (alts); if (_strnicmp (alts, rp, len) == 0) break; alts = strchr (alts, '\0') + 1; ++val; } if (*alts == '\0') { if (*decided == loc && val != 0) return NULL; } else { *decided = loc; break; } } do { val = 0; while (*rp == ' ') ++rp; if (*rp < '0' || *rp > '9') return NULL; do { val *= 10; val += *rp++ - '0'; } while (val * 10 <= to && *rp >= '0' && *rp <= '9'); if (val < from || val > to) return NULL; } while (0); } while (0)
	//#else
	//C++ TO JAVA CONVERTER TODO TASK: Alternate #define macros with the same name cannot be converted to Java:
	//#define get_alt_number(from, to) do { val = 0; while (*rp == ' ') ++rp; if (*rp < '0' || *rp > '9') return NULL; do { val *= 10; val += *rp++ - '0'; } while (val * 10 <= to && *rp >= '0' && *rp <= '9'); if (val < from || val > to) return NULL; } while (0)
	//#endif
	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//#define recursive(new_fmt) (*(new_fmt) != '\0' && (rp = strptime_internal (rp, (new_fmt), tm, decided)) != NULL)


	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if _LIBC
	// This is defined in locale/C-time.c in the GNU libc.  
//C++ TO JAVA CONVERTER NOTE: 'extern' variable declarations are not required in Java:
	//extern const struct locale_data _nl_C_LC_TIME;
//C++ TO JAVA CONVERTER NOTE: 'extern' variable declarations are not required in Java:
	//extern const short __mon_yday[2][13];
//C++ TO JAVA CONVERTER TODO TASK: The following line could not be converted:
static char const (&_nl_C_LC_TIME.values[_NL_ITEM_INDEX (DAY_1)].string)[][10] = { "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday" };
//C++ TO JAVA CONVERTER TODO TASK: The following line could not be converted:
static char const (&_nl_C_LC_TIME.values[_NL_ITEM_INDEX (ABDAY_1)].string)[][4] = { "Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat" };
//C++ TO JAVA CONVERTER TODO TASK: The following line could not be converted:
static char const (&_nl_C_LC_TIME.values[_NL_ITEM_INDEX (MON_1)].string)[][10] = { "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December" };
//C++ TO JAVA CONVERTER TODO TASK: The following line could not be converted:
static char const (&_nl_C_LC_TIME.values[_NL_ITEM_INDEX (ABMON_1)].string)[][4] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };

	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//#define weekday_name (&_nl_C_LC_TIME.values[_NL_ITEM_INDEX (DAY_1)].string)
	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//#define ab_weekday_name (&_nl_C_LC_TIME.values[_NL_ITEM_INDEX (ABDAY_1)].string)
	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//#define month_name (&_nl_C_LC_TIME.values[_NL_ITEM_INDEX (MON_1)].string)
	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//#define ab_month_name (&_nl_C_LC_TIME.values[_NL_ITEM_INDEX (ABMON_1)].string)
	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//#define HERE_D_T_FMT (_nl_C_LC_TIME.values[_NL_ITEM_INDEX (D_T_FMT)].string)
	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//#define HERE_D_FMT (_nl_C_LC_TIME.values[_NL_ITEM_INDEX (D_FMT)].string)
	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//#define HERE_AM_STR (_nl_C_LC_TIME.values[_NL_ITEM_INDEX (AM_STR)].string)
	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//#define HERE_PM_STR (_nl_C_LC_TIME.values[_NL_ITEM_INDEX (PM_STR)].string)
	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//#define HERE_T_FMT_AMPM (_nl_C_LC_TIME.values[_NL_ITEM_INDEX (T_FMT_AMPM)].string)
	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//#define HERE_T_FMT (_nl_C_LC_TIME.values[_NL_ITEM_INDEX (T_FMT)].string)

	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//#define strncasecmp(s1, s2, n) __strncasecmp (s1, s2, n)
	//#else
	//#define HERE_D_T_FMT "%a %b %e %H:%M:%S %Y"
	//#define HERE_D_FMT "%m/%d/%y"
	//#define HERE_AM_STR "AM"
	//#define HERE_PM_STR "PM"
	//#define HERE_T_FMT_AMPM "%I:%M:%S %p"
	//#define HERE_T_FMT "%H:%M:%S"

	// Bert Kampes... = 
	// const unsigned short int __mon_yday[1][13] = 
		// Normal years.  
		// Leap years.  
	public static final short[][] __mon_yday = { { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 }, { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 } };
	//#endif

	// Compute the day of the week.  
	static void day_of_the_week(RefObject<tm> tm)
	{
	//   We know that January 1st 1970 was a Thursday (= 4).  Compute the
	//     the difference between this data in the one on TM and so determine
	//     the weekday.  
	  int corr_year = 1900 + tm.argvalue.tm_year - (tm.argvalue.tm_mon < 2);
	  int wday = (-473 + (365 * (tm.argvalue.tm_year - 70)) + (corr_year / 4) - ((corr_year / 4) / 25) + ((corr_year / 4) % 25 < 0) + (((corr_year / 4) / 25) / 4) + __mon_yday[0][tm.argvalue.tm_mon] + tm.argvalue.tm_mday - 1);
	  tm.argvalue.tm_wday = wday % 7;
	}

	// Compute the day of the year.  
	static void day_of_the_year(RefObject<tm> tm)
	{
	  tm.argvalue.tm_yday = (__mon_yday[((1900 + tm.argvalue.tm_year) % 4 == 0 && ((1900 + tm.argvalue.tm_year) % 100 != 0 || (1900 + tm.argvalue.tm_year) % 400 == 0))][tm.argvalue.tm_mon] + (tm.argvalue.tm_mday - 1));
	}
//C++ TO JAVA CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in Java):
private static char * strptime_internal __P((String buf, String format, struct tm *tm, enum locale_status *decided))


	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if _LIBC
	static char *internal_function strptime_internal __P((String buf, String format, struct tm *tm, enum locale_status *decided))
	//#else
	//C++ TO JAVA CONVERTER NOTE: This static local variable declaration (not allowed in Java) has been moved just prior to the method:
	//static char * strptime_internal __P((String buf, String format, struct tm *tm, enum locale_status *decided))
	//#endif
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro __P was defined in alternate ways and cannot be replaced in-line:
	{
	  public static String rp;
	  public static String fmt;
	  public static int cnt;
	  public static size_t val;
	  public static int have_I;
	  public static int is_pm;
	  public static int century;
	  public static int want_century;
	  public static int have_wday;
	  public static int want_xday;
	  public static int have_yday;

//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	  rp = buf;
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	  fmt = format;
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	  have_I = is_pm = 0;
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	  century = -1;
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	  want_century = 0;
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	  have_wday = want_xday = have_yday = 0;

//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	  while (*fmt != '\0')
		{
	//       A white space in the format string matches 0 more or white
	//         space in the input string.  
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
		  if (isspace (*fmt))
			{
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			  while (isspace (*rp))
				 public static ++rp;
			   public static ++fmt;
			  continue;
			}

	//       Any character but `%' must be matched by the same character
	//         in the iput string.  
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
		  if (*fmt != '%')
			{
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			  if (*fmt++ != *rp++)
				  public static return null;
			  continue;
			}

		   public static ++fmt;
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if ! _NL_CURRENT
		  // We need this for handling the `E' modifier.  
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
		start_over:
	//#endif
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
		  switch (*fmt++)
			{
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			case '%':
			  // Match the `%' character itself.  
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			  if ('%' != *rp++)
				  public static return null;
			  break;
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			case 'a':
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			case 'A':
			  // Match day of week.  
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			  for (cnt = 0; cnt < 7; ++cnt)
				{
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if _NL_CURRENT
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
				  if (*decided !=raw)
					{
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro match_string was defined in alternate ways and cannot be replaced in-line:
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
					  if (match_string(_NL_CURRENT (LC_TIME, DAY_1 + cnt), rp))
						{
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
						  if (*decided == notyet && strcmp (_NL_CURRENT (LC_TIME, DAY_1 + cnt), (&_nl_C_LC_TIME.values[_NL_ITEM_INDEX (DAY_1)].string)[cnt]))
							public static * decided = locale_status.loc;
						  break;
						}
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro match_string was defined in alternate ways and cannot be replaced in-line:
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
					  if (match_string(_NL_CURRENT (LC_TIME, ABDAY_1 + cnt), rp))
						{
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
						  if (*decided == notyet && strcmp (_NL_CURRENT (LC_TIME, ABDAY_1 + cnt), (&_nl_C_LC_TIME.values[_NL_ITEM_INDEX (ABDAY_1)].string)[cnt]))
							public static * decided = locale_status.loc;
						  break;
						}
					}
	//#endif
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro match_string was defined in alternate ways and cannot be replaced in-line:
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
				  if (*decided != loc && (match_string((&_nl_C_LC_TIME.values[_NL_ITEM_INDEX (DAY_1)].string)[cnt], rp) || match_string((&_nl_C_LC_TIME.values[_NL_ITEM_INDEX (ABDAY_1)].string)[cnt], rp)))
					{
					  public static * decided = locale_status.raw;
					  break;
					}
				}
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			  if (cnt == 7)
				// Does notyet match a weekday name.  
				public static return null;
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			  tm.tm_wday = cnt;
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			  have_wday = 1;
			  break;
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			case 'b':
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			case 'B':
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			case 'h':
			  // Match month name.  
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			  for (cnt = 0; cnt < 12; ++cnt)
				{
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if _NL_CURRENT
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
				  if (*decided !=raw)
					{
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro match_string was defined in alternate ways and cannot be replaced in-line:
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
					  if (match_string(_NL_CURRENT (LC_TIME, MON_1 + cnt), rp))
						{
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
						  if (*decided == notyet && strcmp (_NL_CURRENT (LC_TIME, MON_1 + cnt), (&_nl_C_LC_TIME.values[_NL_ITEM_INDEX (MON_1)].string)[cnt]))
							public static * decided = locale_status.loc;
						  break;
						}
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro match_string was defined in alternate ways and cannot be replaced in-line:
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
					  if (match_string(_NL_CURRENT (LC_TIME, ABMON_1 + cnt), rp))
						{
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
						  if (*decided == notyet && strcmp (_NL_CURRENT (LC_TIME, ABMON_1 + cnt), (&_nl_C_LC_TIME.values[_NL_ITEM_INDEX (ABMON_1)].string)[cnt]))
							public static * decided = locale_status.loc;
						  break;
						}
					}
	//#endif
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro match_string was defined in alternate ways and cannot be replaced in-line:
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
				  if (match_string((&_nl_C_LC_TIME.values[_NL_ITEM_INDEX (MON_1)].string)[cnt], rp) || match_string((&_nl_C_LC_TIME.values[_NL_ITEM_INDEX (ABMON_1)].string)[cnt], rp))
					{
					  public static * decided = locale_status.raw;
					  break;
					}
				}
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			  if (cnt == 12)
				// Does notyet match a month name.  
				public static return null;
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			  tm.tm_mon = cnt;
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			  want_xday = 1;
			  break;
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			case 'c':
			  // Match locale's date and time format.  
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if _NL_CURRENT
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
			  if (*decided != raw)
				{
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
				  if (!(*(_NL_CURRENT (LC_TIME, D_T_FMT)) != '\0' && (rp = strptime_internal (rp, (_NL_CURRENT (LC_TIME, D_T_FMT)), tm, decided)) != null))
					{
//C++ TO JAVA CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
					  if (*decided == loc)
						public static return null;
					}
				  else
					{
					  if (*decided == locale_status.notyet && strcmp (_NL_CURRENT (LC_TIME, D_T_FMT), (_nl_C_LC_TIME.values[_NL_ITEM_INDEX (D_T_FMT)].string)))
						*decided = locale_status.loc;
					  want_xday = 1;
					  break;
					}
				  *decided = locale_status.raw;
				}
	//#endif
			  if (!(*((_nl_C_LC_TIME.values[_NL_ITEM_INDEX (D_T_FMT)].string)) != '\0' && (rp = strptime_internal (rp, ((_nl_C_LC_TIME.values[_NL_ITEM_INDEX (D_T_FMT)].string)), tm, decided)) != null))
				return null;
			  want_xday = 1;
			  break;
			case 'C':
			  // Match century number.  
			  do
			  {
				  val = 0;
				  while (*rp == ' ')
					  ++rp;
				  if (*rp < '0' || *rp > '9')
					  return null;
				  do
				  {
					  val *= 10;
					  val += *rp++ - '0';
				  } while (val * 10 <= 99 && *rp >= '0' && *rp <= '9');
				  if (val < 0 || val > 99)
					  return null;
				  } while (0);
			  century = val;
			  want_xday = 1;
			  break;
			case 'd':
			case 'e':
			  // Match day of month.  
			  do
			  {
				  val = 0;
				  while (*rp == ' ')
					  ++rp;
				  if (*rp < '0' || *rp > '9')
					  return null;
				  do
				  {
					  val *= 10;
					  val += *rp++ - '0';
				  } while (val * 10 <= 31 && *rp >= '0' && *rp <= '9');
				  if (val < 1 || val > 31)
					  return null;
				  } while (0);
			  tm.tm_mday = val;
			  want_xday = 1;
			  break;
			case 'F':
			  if (!(*("%Y-%m-%d") != '\0' && (rp = strptime_internal (rp, ("%Y-%m-%d"), tm, decided)) != null))
				return null;
			  want_xday = 1;
			  break;
			case 'x':
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if _NL_CURRENT
			  if (*decided != locale_status.raw)
				{
				  if (!(*(_NL_CURRENT (LC_TIME, D_FMT)) != '\0' && (rp = strptime_internal (rp, (_NL_CURRENT (LC_TIME, D_FMT)), tm, decided)) != null))
					{
					  if (*decided == locale_status.loc)
						return null;
					}
				  else
					{
					  if (decided == locale_status.notyet && strcmp (_NL_CURRENT (LC_TIME, D_FMT), (_nl_C_LC_TIME.values[_NL_ITEM_INDEX (D_FMT)].string)))
						*decided = locale_status.loc;
					  want_xday = 1;
					  break;
					}
				  *decided = locale_status.raw;
				}
	//#endif
			  // Fall through.  
			case 'D':
			  // Match standard day format.  
			  if (!(*((_nl_C_LC_TIME.values[_NL_ITEM_INDEX (D_FMT)].string)) != '\0' && (rp = strptime_internal (rp, ((_nl_C_LC_TIME.values[_NL_ITEM_INDEX (D_FMT)].string)), tm, decided)) != null))
				return null;
			  want_xday = 1;
			  break;
			case 'k':
			case 'H':
			  // Match hour in 24-hour clock.  
			  do
			  {
				  val = 0;
				  while (*rp == ' ')
					  ++rp;
				  if (*rp < '0' || *rp > '9')
					  return null;
				  do
				  {
					  val *= 10;
					  val += *rp++ - '0';
				  } while (val * 10 <= 23 && *rp >= '0' && *rp <= '9');
				  if (val < 0 || val > 23)
					  return null;
				  } while (0);
			  tm.tm_hour = val;
			  have_I = 0;
			  break;
			case 'I':
			  // Match hour in 12-hour clock.  
			  do
			  {
				  val = 0;
				  while (*rp == ' ')
					  ++rp;
				  if (*rp < '0' || *rp > '9')
					  return null;
				  do
				  {
					  val *= 10;
					  val += *rp++ - '0';
				  } while (val * 10 <= 12 && *rp >= '0' && *rp <= '9');
				  if (val < 1 || val > 12)
					  return null;
				  } while (0);
			  tm.tm_hour = val % 12;
			  have_I = 1;
			  break;
			case 'j':
			  // Match day number of year.  
			  do
			  {
				  val = 0;
				  while (*rp == ' ')
					  ++rp;
				  if (*rp < '0' || *rp > '9')
					  return null;
				  do
				  {
					  val *= 10;
					  val += *rp++ - '0';
				  } while (val * 10 <= 366 && *rp >= '0' && *rp <= '9');
				  if (val < 1 || val > 366)
					  return null;
				  } while (0);
			  tm.tm_yday = val - 1;
			  have_yday = 1;
			  break;
			case 'm':
			  // Match number of month.  
			  do
			  {
				  val = 0;
				  while (*rp == ' ')
					  ++rp;
				  if (*rp < '0' || *rp > '9')
					  return null;
				  do
				  {
					  val *= 10;
					  val += *rp++ - '0';
				  } while (val * 10 <= 12 && *rp >= '0' && *rp <= '9');
				  if (val < 1 || val > 12)
					  return null;
				  } while (0);
			  tm.tm_mon = val - 1;
			  want_xday = 1;
			  break;
			case 'M':
			  // Match minute.  
			  do
			  {
				  val = 0;
				  while (*rp == ' ')
					  ++rp;
				  if (*rp < '0' || *rp > '9')
					  return null;
				  do
				  {
					  val *= 10;
					  val += *rp++ - '0';
				  } while (val * 10 <= 59 && *rp >= '0' && *rp <= '9');
				  if (val < 0 || val > 59)
					  return null;
				  } while (0);
			  tm.tm_min = val;
			  break;
			case 'n':
			case 't':
			  // Match any white space.  
			  while (Character.isWhitespace (*rp))
				++rp;
			  break;
			case 'p':
			  // Match locale's equivalent of AM/PM.  
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if _NL_CURRENT
			  if (*decided != locale_status.raw)
				{
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro match_string was defined in alternate ways and cannot be replaced in-line:
				  if (match_string(_NL_CURRENT (LC_TIME, AM_STR), rp))
					{
					  if (strcmp (_NL_CURRENT (LC_TIME, AM_STR), (_nl_C_LC_TIME.values[_NL_ITEM_INDEX (AM_STR)].string)))
						*decided = locale_status.loc;
					  break;
					}
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro match_string was defined in alternate ways and cannot be replaced in-line:
				  if (match_string(_NL_CURRENT (LC_TIME, PM_STR), rp))
					{
					  if (strcmp (_NL_CURRENT (LC_TIME, PM_STR), (_nl_C_LC_TIME.values[_NL_ITEM_INDEX (PM_STR)].string)))
						*decided = locale_status.loc;
					  is_pm = 1;
					  break;
					}
				  *decided = locale_status.raw;
				}
	//#endif
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro match_string was defined in alternate ways and cannot be replaced in-line:
			  if (!match_string((_nl_C_LC_TIME.values[_NL_ITEM_INDEX (AM_STR)].string), rp))
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro match_string was defined in alternate ways and cannot be replaced in-line:
				if (match_string((_nl_C_LC_TIME.values[_NL_ITEM_INDEX (PM_STR)].string), rp))
				  is_pm = 1;
				else
				  return null;
			  break;
			case 'r':
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if _NL_CURRENT
			  if (*decided != locale_status.raw)
				{
				  if (!(*(_NL_CURRENT (LC_TIME, T_FMT_AMPM)) != '\0' && (rp = strptime_internal (rp, (_NL_CURRENT (LC_TIME, T_FMT_AMPM)), tm, decided)) != null))
					{
					  if (*decided == locale_status.loc)
						return null;
					}
				  else
					{
					  if (*decided == locale_status.notyet && strcmp (_NL_CURRENT (LC_TIME, T_FMT_AMPM), (_nl_C_LC_TIME.values[_NL_ITEM_INDEX (T_FMT_AMPM)].string)))
						*decided = locale_status.loc;
					  break;
					}
				  *decided = locale_status.raw;
				}
	//#endif
			  if (!(*((_nl_C_LC_TIME.values[_NL_ITEM_INDEX (T_FMT_AMPM)].string)) != '\0' && (rp = strptime_internal (rp, ((_nl_C_LC_TIME.values[_NL_ITEM_INDEX (T_FMT_AMPM)].string)), tm, decided)) != null))
				return null;
			  break;
			case 'R':
			  if (!(*("%H:%M") != '\0' && (rp = strptime_internal (rp, ("%H:%M"), tm, decided)) != null))
				return null;
			  break;
			case 's':
			  {
	//             The number of seconds may be very high so we cannot use
	//               the `get_number' macro.  Instead read the number
	//               character for character and construct the result while
	//               doing this.  
				time_t secs = 0;
				if (*rp < '0' || *rp > '9')
				  // We need at least one digit.  
				  return null;

				do
				  {
					secs *= 10;
					secs += *rp++ - '0';
				  } while (*rp >= '0' && *rp <= '9');

	//C++ TO JAVA CONVERTER TODO TASK: The #define macro localtime_r was defined in alternate ways and cannot be replaced in-line:
				if (localtime_r (secs, tm) == null)
				  // Error in function.  
				  return null;
			  }
			  break;
			case 'S':
			  do
			  {
				  val = 0;
				  while (*rp == ' ')
					  ++rp;
				  if (*rp < '0' || *rp > '9')
					  return null;
				  do
				  {
					  val *= 10;
					  val += *rp++ - '0';
				  } while (val * 10 <= 61 && *rp >= '0' && *rp <= '9');
				  if (val < 0 || val > 61)
					  return null;
				  } while (0);
			  tm.tm_sec = val;
			  break;
			case 'X':
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if _NL_CURRENT
			  if (*decided != locale_status.raw)
				{
				  if (!(*(_NL_CURRENT (LC_TIME, T_FMT)) != '\0' && (rp = strptime_internal (rp, (_NL_CURRENT (LC_TIME, T_FMT)), tm, decided)) != null))
					{
					  if (*decided == locale_status.loc)
						return null;
					}
				  else
					{
					  if (strcmp (_NL_CURRENT (LC_TIME, T_FMT), (_nl_C_LC_TIME.values[_NL_ITEM_INDEX (T_FMT)].string)))
						*decided = locale_status.loc;
					  break;
					}
				  *decided = locale_status.raw;
				}
	//#endif
			  // Fall through.  
			case 'T':
			  if (!(*((_nl_C_LC_TIME.values[_NL_ITEM_INDEX (T_FMT)].string)) != '\0' && (rp = strptime_internal (rp, ((_nl_C_LC_TIME.values[_NL_ITEM_INDEX (T_FMT)].string)), tm, decided)) != null))
				return null;
			  break;
			case 'u':
			  do
			  {
				  val = 0;
				  while (*rp == ' ')
					  ++rp;
				  if (*rp < '0' || *rp > '9')
					  return null;
				  do
				  {
					  val *= 10;
					  val += *rp++ - '0';
				  } while (val * 10 <= 7 && *rp >= '0' && *rp <= '9');
				  if (val < 1 || val > 7)
					  return null;
				  } while (0);
			  tm.tm_wday = val % 7;
			  have_wday = 1;
			  break;
			case 'g':
			  do
			  {
				  val = 0;
				  while (*rp == ' ')
					  ++rp;
				  if (*rp < '0' || *rp > '9')
					  return null;
				  do
				  {
					  val *= 10;
					  val += *rp++ - '0';
				  } while (val * 10 <= 99 && *rp >= '0' && *rp <= '9');
				  if (val < 0 || val > 99)
					  return null;
				  } while (0);
			  // XXX This cannot determine any field in TM.  
			  break;
			case 'G':
			  if (*rp < '0' || *rp > '9')
				return null;
	//           XXX Ignore the number since we would need some more
	//             information to compute a real date.  
			  do
			  {
				++rp;
			  }while (*rp >= '0' && *rp <= '9');
			  break;
			case 'U':
			case 'V':
			case 'W':
			  do
			  {
				  val = 0;
				  while (*rp == ' ')
					  ++rp;
				  if (*rp < '0' || *rp > '9')
					  return null;
				  do
				  {
					  val *= 10;
					  val += *rp++ - '0';
				  } while (val * 10 <= 53 && *rp >= '0' && *rp <= '9');
				  if (val < 0 || val > 53)
					  return null;
				  } while (0);
	//           XXX This cannot determine any field in TM without some
	//             information.  
			  break;
			case 'w':
			  // Match number of weekday.  
			  do
			  {
				  val = 0;
				  while (*rp == ' ')
					  ++rp;
				  if (*rp < '0' || *rp > '9')
					  return null;
				  do
				  {
					  val *= 10;
					  val += *rp++ - '0';
				  } while (val * 10 <= 6 && *rp >= '0' && *rp <= '9');
				  if (val < 0 || val > 6)
					  return null;
				  } while (0);
			  tm.tm_wday = val;
			  have_wday = 1;
			  break;
			case 'y':
			  // Match year within century.  
			  do
			  {
				  val = 0;
				  while (*rp == ' ')
					  ++rp;
				  if (*rp < '0' || *rp > '9')
					  return null;
				  do
				  {
					  val *= 10;
					  val += *rp++ - '0';
				  } while (val * 10 <= 99 && *rp >= '0' && *rp <= '9');
				  if (val < 0 || val > 99)
					  return null;
				  } while (0);
	//           The "Year 2000: The Millennium Rollover" paper suggests that
	//             values in the range 69-99 refer to the twentieth century.  
			  tm.tm_year = val >= 69 ? val : val + 100;
			  // Indicate that we want to use the century, if specified.  
			  want_century = 1;
			  want_xday = 1;
			  break;
			case 'Y':
			  // Match year including century number.  
			  do
			  {
				  val = 0;
				  while (*rp == ' ')
					  ++rp;
				  if (*rp < '0' || *rp > '9')
					  return null;
				  do
				  {
					  val *= 10;
					  val += *rp++ - '0';
				  } while (val * 10 <= 9999 && *rp >= '0' && *rp <= '9');
				  if (val < 0 || val > 9999)
					  return null;
				  } while (0);
			  tm.tm_year = val - 1900;
			  want_century = 0;
			  want_xday = 1;
			  break;
			case 'Z':
			  // XXX How to handle this?  
			  break;
			case 'E':
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if _NL_CURRENT
			  switch (*fmt++)
				{
				case 'c':
				  // Match locale's alternate date and time format.  
				  if (*decided != locale_status.raw)
					{
					  String fmt = _NL_CURRENT (LC_TIME, ERA_D_T_FMT);

					  if ( fmt.equals('\0'))
						fmt = _NL_CURRENT (LC_TIME, D_T_FMT);

					  if (!(*(fmt) != '\0' && (rp = strptime_internal (rp, (fmt), tm, decided)) != null))
						{
						  if (*decided == locale_status.loc)
							return null;
						}
					  else
						{
						  if (strcmp (fmt, (_nl_C_LC_TIME.values[_NL_ITEM_INDEX (D_T_FMT)].string)))
							*decided = locale_status.loc;
						  want_xday = 1;
						  break;
						}
					  *decided = locale_status.raw;
					}
	//               The C locale has no era information, so use the
	//                 normal representation.  
				  if (!(*((_nl_C_LC_TIME.values[_NL_ITEM_INDEX (D_T_FMT)].string)) != '\0' && (rp = strptime_internal (rp, ((_nl_C_LC_TIME.values[_NL_ITEM_INDEX (D_T_FMT)].string)), tm, decided)) != null))
					return null;
				  want_xday = 1;
				  break;
				case 'C':
				case 'y':
				case 'Y':
	//               Match name of base year in locale's alternate
	//                 representation.  
	//               XXX This is currently not implemented.  It should
	//                 use the value _NL_CURRENT (LC_TIME, ERA).  
				  break;
				case 'x':
				  if (*decided != locale_status.raw)
					{
					  String fmt = _NL_CURRENT (LC_TIME, ERA_D_FMT);

					  if ( fmt.equals('\0'))
						fmt = _NL_CURRENT (LC_TIME, D_FMT);

					  if (!(*(fmt) != '\0' && (rp = strptime_internal (rp, (fmt), tm, decided)) != null))
						{
						  if (*decided == locale_status.loc)
							return null;
						}
					  else
						{
						  if (strcmp (fmt, (_nl_C_LC_TIME.values[_NL_ITEM_INDEX (D_FMT)].string)))
							*decided = locale_status.loc;
						  break;
						}
					  *decided = locale_status.raw;
					}
				  if (!(*((_nl_C_LC_TIME.values[_NL_ITEM_INDEX (D_FMT)].string)) != '\0' && (rp = strptime_internal (rp, ((_nl_C_LC_TIME.values[_NL_ITEM_INDEX (D_FMT)].string)), tm, decided)) != null))
					return null;
				  break;
				case 'X':
				  if (*decided != locale_status.raw)
					{
					  String fmt = _NL_CURRENT (LC_TIME, ERA_T_FMT);

					  if ( fmt.equals('\0'))
						fmt = _NL_CURRENT (LC_TIME, T_FMT);

					  if (!(*(fmt) != '\0' && (rp = strptime_internal (rp, (fmt), tm, decided)) != null))
						{
						  if (*decided == locale_status.loc)
							return null;
						}
					  else
						{
						  if (strcmp (fmt, (_nl_C_LC_TIME.values[_NL_ITEM_INDEX (T_FMT)].string)))
							*decided = locale_status.loc;
						  break;
						}
					  *decided = locale_status.raw;
					}
				  if (!(*((_nl_C_LC_TIME.values[_NL_ITEM_INDEX (T_FMT)].string)) != '\0' && (rp = strptime_internal (rp, ((_nl_C_LC_TIME.values[_NL_ITEM_INDEX (T_FMT)].string)), tm, decided)) != null))
					return null;
				  break;
				default:
				  return null;
				}
			  break;
	//#else
	//           We have no information about the era format.  Just use
	//             the normal format.  
			  if ( !fmt.equals('c') && !fmt.equals('C') && !fmt.equals('y') && !fmt.equals('Y') && !fmt.equals('x') && !fmt.equals('X'))
				// This is an illegal format.  
				return null;

//C++ TO JAVA CONVERTER TODO TASK: There are no gotos or labels in Java:
			  goto start_over;
	//#endif
			case 'O':
			  switch ( fmt++)
				{
				case 'd':
				case 'e':
				  // Match day of month using alternate numeric symbols.  
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro get_alt_number was defined in alternate ways and cannot be replaced in-line:
				  get_alt_number(1, 31);
				  tm.tm_mday = val;
				  want_xday = 1;
				  break;
				case 'H':
	//               Match hour in 24-hour clock using alternate numeric
	//                 symbols.  
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro get_alt_number was defined in alternate ways and cannot be replaced in-line:
				  get_alt_number(0, 23);
				  tm.tm_hour = val;
				  have_I = 0;
				  break;
				case 'I':
	//               Match hour in 12-hour clock using alternate numeric
	//                 symbols.  
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro get_alt_number was defined in alternate ways and cannot be replaced in-line:
				  get_alt_number(1, 12);
				  tm.tm_hour = val - 1;
				  have_I = 1;
				  break;
				case 'm':
				  // Match month using alternate numeric symbols.  
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro get_alt_number was defined in alternate ways and cannot be replaced in-line:
				  get_alt_number(1, 12);
				  tm.tm_mon = val - 1;
				  want_xday = 1;
				  break;
				case 'M':
				  // Match minutes using alternate numeric symbols.  
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro get_alt_number was defined in alternate ways and cannot be replaced in-line:
				  get_alt_number(0, 59);
				  tm.tm_min = val;
				  break;
				case 'S':
				  // Match seconds using alternate numeric symbols.  
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro get_alt_number was defined in alternate ways and cannot be replaced in-line:
				  get_alt_number(0, 61);
				  tm.tm_sec = val;
				  break;
				case 'U':
				case 'V':
				case 'W':
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro get_alt_number was defined in alternate ways and cannot be replaced in-line:
				  get_alt_number(0, 53);
	//               XXX This cannot determine any field in TM without
	//                 further information.  
				  break;
				case 'w':
				  // Match number of weekday using alternate numeric symbols.  
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro get_alt_number was defined in alternate ways and cannot be replaced in-line:
				  get_alt_number(0, 6);
				  tm.tm_wday = val;
				  have_wday = 1;
				  break;
				case 'y':
				  // Match year within century using alternate numeric symbols.  
	//C++ TO JAVA CONVERTER TODO TASK: The #define macro get_alt_number was defined in alternate ways and cannot be replaced in-line:
				  get_alt_number(0, 99);
				  tm.tm_year = val >= 69 ? val : val + 100;
				  want_xday = 1;
				  break;
				default:
				  return null;
				}
			  break;
			default:
			  return null;
			}
		}

	  if (have_I && is_pm)
		tm.tm_hour += 12;

	  if (want_century && century != -1)
		tm.tm_year = tm.tm_year % 100 + (century - 19) * 100;

	  if (want_xday && !have_wday)
		day_of_the_week (tm);
	  if (want_xday && !have_yday)
		day_of_the_year (tm);

	  return (char) rp;
	}


	char * strptime (String buf, String format, struct tm *tm)
	{
	  locale_status decided;
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if _NL_CURRENT
	  decided = locale_status.notyet;
	//#else
	  decided = locale_status.raw;
	//#endif
	  return strptime_internal (buf, format, tm, decided);
	}


	//#endif // __NO_STRPTIME

}
//
// * Copyright (c) 1999-2009 Delft University of Technology, The Netherlands
// *
// * This file is part of Doris, the Delft o-o radar interferometric software.
// *
// * Doris program is free software; you can redistribute it and/or modify
// * it under the terms of the GNU General Public License as published by
// * the Free Software Foundation; either version 2 of the License, or
// * (at your option) any later version.
// *
// * Doris is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// * GNU General Public License for more details.
// *
// * You should have received a copy of the GNU General Public License
// * along with this program; if not, write to the Free Software
// * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// *
// * Publications that contain results produced by the Doris software should
// * contain an acknowledgment. (For example: The interferometric processing
// * was performed using the freely available Doris software package developed
// * by the Delft Institute of Earth Observation and Space Systems (DEOS), Delft
// * University of Technology, or include a reference to: Bert Kampes and
// * Stefania Usai. \"Doris: The Delft Object-oriented Radar Interferometric
// * software.\" In: proceedings 2nd ITC ORS symposium, August 1999. (cdrom)).
// *
// 
//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
//#if __NO_STRPTIME
//
//The function "strptime" was unknown on my linux system.
//I downloaded the source, made small changes in order to compile
//it with g++, and include it with the Doris software.
//I did not perform a good check.
//DO NOT USE THIS, unless you cannot compile the doris software
//otherwise. in that case, define __NO_STRPTIME, in the makefile.
//That will use this code.
//Bert Kampes, November 2000.
//
//
//time/strptime.c
//compile: g++ strptime.cc teststrptime.c
//
//FUNCTIONS
//This source file includes following functions. 
//localtime_r 
//match_char 
//match_string 
//match_string 
//get_number 
//get_alt_number 
//get_alt_number 
//recursive 
//strncasecmp 
//__isleap 
//day_of_the_week 
//day_of_the_year 
//strptime_internal 
//strptime 
//

// Convert a string representation of time to a time value.
//   Copyright (C) 1996, 1997, 1998, 1999 Free Software Foundation, Inc.
//   This file is part of the GNU C Library.
//   Contributed by Ulrich Drepper <drepper@cygnus.com>, 1996.
//
//   The GNU C Library is free software; you can redistribute it and/or
//   modify it under the terms of the GNU Library General Public License as
//   published by the Free Software Foundation; either version 2 of the
//   License, or (at your option) any later version.
//
//   The GNU C Library is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//   Library General Public License for more details.
//
//   You should have received a copy of the GNU Library General Public
//   License along with the GNU C Library; see the file COPYING.LIB.  If not,
//   write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
//   Boston, MA 02111-1307, USA.  

// XXX This version of the implementation is not really complete.
//   Some of the fields cannot add information alone.  But if seeing
//   some of them in the same format (such as year, week and weekday)
//   this is enough information for determining the date.  

// #ifdef HAVE_CONFIG_H
// #include <config.h>
// #endif

//#include <ctype.h>
// BK 03-Nov-2002: changed include files to cctype etc. not ctype.h
// commented this out, BK, seems not required and fails under cygwin.
// #include <langinfo.h>

//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
//#if _LIBC
//# include "../locale/localeinfo.h"
//#endif


//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
//#if ! __P
//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
//#if __GNUC__ || (__STDC__ && __STDC__)
//C++ TO JAVA CONVERTER TODO TASK: Alternate #define macros with the same name cannot be converted to Java:
//#define __P(args) args
//#else
//C++ TO JAVA CONVERTER TODO TASK: Alternate #define macros with the same name cannot be converted to Java:
//#define __P(args) ()
//#endif // GCC.
//#endif // Not __P.

// Status of lookup: do we use the locale data or the raw data?  
public enum locale_status
{
	notyet,
	loc,
	raw
}


//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
//#if ! __isleap
// Nonzero if YEAR is a leap year (every 4 years,
//   except every 100th isn't, and every 400th is).  
//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
//#define __isleap(year) ((year) % 4 == 0 && ((year) % 100 != 0 || (year) % 400 == 0))
//#endif

final class DefineConstantsTmp_strptime
{
	public static final String HERE_D_T_FMT = "%a %b %e %H:%M:%S %Y";
	public static final String HERE_D_FMT = "%m/%d/%y";
	public static final String HERE_AM_STR = "AM";
	public static final String HERE_PM_STR = "PM";
	public static final String HERE_T_FMT_AMPM = "%I:%M:%S %p";
	public static final String HERE_T_FMT = "%H:%M:%S";
}
//----------------------------------------------------------------------------------------
//	Copyright © 2006 - 2008 Tangible Software Solutions Inc.
//
//	This class is used to simulate the ability to pass arguments by reference in Java.
//----------------------------------------------------------------------------------------
final class RefObject<T>
{
	T argvalue;
	RefObject(T refarg)
	{
		argvalue = refarg;
	}
}