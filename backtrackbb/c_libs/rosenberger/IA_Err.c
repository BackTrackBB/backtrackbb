/*------------------------------------------------------------*
 * IA_Log_Die / IA_Log_Warn
 *
 * Log to Syslog, warn or log and commit suicide
 *
 * ALso provides seismic event logging to syslog user.alert
 * in IA_Log_Event
 *------------------------------------------------------------*/

#include <stdio.h>
#include <unistd.h>
#include <syslog.h>
#include <stdarg.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include "IA_Kdiag.h"
#include "IA_Err.h"

#define IA_ERRFILE "PS_FiltErr.log"
#define IA_DATAPATH "./"

/*------------------------------------------------------------*/

/** @name IA_Log_Die
 * Writes error message to syslog LOG_DAMON,
 * with priority LOG_ALERT.
 * Signals all threads to cleanup
 * and terminate, calls exit(1)
 * @param mesg : Error message to log, variable parameter list
 * @return: It does not return
 */
void IA_Log_Die(char *mesg, ...)
{
    va_list ap;
    FILE *ep;
    char *fn = IA_ERRFILE;

#ifndef IA_DEBUG
    char buf[1024];
#endif

    va_start(ap,mesg);

    if(strlen(mesg) < 1024) {

        /* try making the directory, case it does not exist yet */
        mkdir(IA_DATAPATH,00777);
        ep = fopen(fn,"w");

#ifndef IA_DEBUG
        vsprintf(buf,mesg,ap);
        openlog(VERSION,LOG_CONS | LOG_NDELAY | LOG_PID,
                LOG_DAEMON );
        syslog(LOG_ALERT,buf);
#else
        vfprintf(stderr,mesg,ap);
        fprintf(stderr,"\n");
#endif
        /* try to write directly into the data directory
         * to alert evrybody */

        if(ep != NULL) {
            vfprintf(ep,mesg,ap);
            fclose(ep);
        }
    }
    va_end(ap);

    exit(1);
}

/*------------------------------------------------------------------*/

/** @name IA_Log_Warn
 * Writes error message to syslog LOG_DAMON,
 * with priority LOG_ALERT.
 * @param mesg : Error message to log, variable parameter list
 * @return nothing
 */
void IA_Log_Warn(const char *mesg, ...)
{
    va_list ap;

#ifndef IA_DEBUG
    char buf[1024];
#endif

    va_start(ap,mesg);

    if(strlen(mesg) < 1024) {
#ifndef IA_DEBUG
        vsprintf(buf,mesg,ap);
        openlog(VERSION,LOG_CONS | LOG_NDELAY | LOG_PID,
                LOG_DAEMON );
        syslog(LOG_NOTICE,buf);
#else
        vfprintf(stderr,mesg,ap);
        fprintf(stderr,"\n");
#endif
    }
    va_end(ap);
}
