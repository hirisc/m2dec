#include <smachine.h>
#include <string.h>
#include "sbrk.h"

#pragma section BOOT

//#define	P_PHI	20000000	/* PÉ”=20MHz */
#define	P_PHI	27000000	/* PÉ”=27MHz */
//#define	P_PHI	33333333	/* PÉ”=33.33MHz */

static void init_CPU(void);							/* CPUÇÃèâä˙âª */
static void LED (int p);
static int tokenize(char *args, char **arg_dst);

extern void init_sbrk();
int main(int argc, char *argv[]);

#define MAX_ARG 16

/***********************************************************************/

/* mainä÷êî */
#pragma entry PowerON_Reset(sp=0x0e000000)

void PowerON_Reset(void)
{
	char args[] = "h264dec dummy.264";
	char *arg[MAX_ARG];
	int arg_num;

	/* CPUèâä˙âª */
	*(volatile long *)0xff00001c = 0x00000105;
	set_cr((get_cr() & ~0x10000000) | 0x1000);
//	set_vbr((unsigned char *)__sectop("PEXPT") - 0x100);			/*   CS0  SRAMãÛä‘ */

	set_cr(get_cr() | 0x1000);
//	_INITSCT();
//	_INIT_IOLIB();
	init_sbrk();	// initialize pointer to the heap area

	arg_num = tokenize(args, arg);
	main(arg_num, &arg);
	sleep();
}


static void LED (int p) {
#ifndef SOFT_ONLY
	*(unsigned short *)0xb0800000 = p;
#endif
}

static int tokenize(char *args, char **arg_dst)
{
	int i;
	char *p;
	char *p_prev = 0;

	for (i = 0, p = strtok(args, " "); p && (i < MAX_ARG); ++i, p = strtok(0, " ")) {
		arg_dst[i] = p;
		if (p_prev) {
			p_prev = '\0';
		}
		p_prev = p;
	}
	return i;
}
