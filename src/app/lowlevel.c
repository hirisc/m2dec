#include <stdlib.h>
#include <stdio.h>

#ifdef __RENESAS_VERSION__

#include <machine.h>

#pragma section FILES

#define O_RDONLY 1
#define O_WRONLY 2
#define MAX_FILESIZE 8 * 1024 * 1024
#define MAX_MD5SIZE 128 * 1024
char InputFile[MAX_FILESIZE];
char OutputFile[MAX_MD5SIZE];
#pragma section
const int InputSize = sizeof(InputFile);
void abort() {while (1);}
static int infilepos;

void exit(int err) {
	while (err) sleep();
}

int fseek(FILE *fp, long offset, int whence) {
	return 0;
}

int open(char *name, int mode, int flg) {
	if (mode & O_RDONLY) {
		infilepos = 0;
	} else if (mode & O_WRONLY) {
		infilepos = 0;
	}
	return 4;
}
int close(int fineno) {return 0;}
int read(int fileno, char *buf, unsigned int count) {
	if (sizeof(InputFile) < infilepos + count) {
		count = sizeof(InputFile) - infilepos;
		if ((signed)count <= 0) {
			return 0;
		}
	}
	memcpy(buf, InputFile + infilepos, count);
	infilepos += count;
	return count;
}
int lseek() {return 0;}
int write(int fileno, char *buf, unsigned int count) {return count;}
char *getenv(const char *name) {return 0;}

#endif /* __RENESAS_VERSION__ */
