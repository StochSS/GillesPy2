
/* Platform-specific macros and functions between POSIX and Windows.
 * Necessary to account for things like interrupt handlers, etc.
 * 
 * Written by: Joshua Cooper
 * April 19, 2021
 */

#ifndef __PLATFORM_H
#define __PLATFORM_H

#ifdef _WIN32
#include <windows.h>

// Interrupt return type, signal number, and return value all differ between Windows and POSIX.
#define INTERRUPT_HANDLER BOOL WINAPI
#define INTERRUPT_SIGNUM DWORD
// Return TRUE on Windows means "this is the last event handler; do not continue."
// This is necessary to ensure that the simulation doesn't stop, so we can cleanly finish simulation.
#define INTERRUPT_RETURN return TRUE

#define SET_INTERRUPT_HANDLER(fn) SetConsoleCtrlHandler(fn, TRUE)

#else
#include <signal.h>

#define INTERRUPT_HANDLER void
#define INTERRUPT_SIGNUM int
#define INTERRUPT_RETURN return

#define SET_INTERRUPT_HANDLER(fn) signal(SIGINT, fn)

#endif

#endif
