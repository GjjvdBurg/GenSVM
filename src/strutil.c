#include "strutil.h"

bool str_startswith(const char *str, const char *pre)
{
	size_t lenpre = strlen(pre),
	       lenstr = strlen(str);
	return lenstr < lenpre ? false : strncmp(pre, str, lenpre) == 0;
}

bool str_endswith(const char *str, const char *suf)
{
	size_t lensuf = strlen(suf),
	       lenstr = strlen(str);
	return lenstr < lensuf ? false : strncmp(str + lenstr - lensuf, suf, lensuf) == 0;
}

void next_line(FILE *fid, char *filename)
{
	char buffer[MAX_LINE_LENGTH];
	get_line(fid, filename, buffer);
}

void get_line(FILE *fid, char *filename, char *buffer)
{
	if (fgets(buffer, MAX_LINE_LENGTH, fid) == NULL) {
		fprintf(stderr, "Error reading file %s\n", filename);
		exit(1);
	}
}

double get_fmt_double(FILE *fid, char *filename, const char *fmt)
{
	char buffer[MAX_LINE_LENGTH];
	double value;

	get_line(fid, filename, buffer);
	sscanf(buffer, fmt, &value);
	return value;
}

long get_fmt_long(FILE *fid, char *filename, const char *fmt)
{
	char buffer[MAX_LINE_LENGTH];
	long value;

	get_line(fid, filename, buffer);
	sscanf(buffer, fmt, &value);
	return value;
}

long all_doubles_str(char *buffer, long offset, double *all_doubles)
{
	double value;
	long i = 0;
	char *start, *end;

	start = buffer + offset;
	while (true) {
		value = strtod(start, &end);
		if (start != end) {
			all_doubles[i] = value;
			i++;
		} else 
			break;
		start = end;
		end = NULL;
	}

	return i;
}

long all_longs_str(char *buffer, long offset, long *all_longs)
{
	long value;
	long i = 0;
	char *start, *end;

	start = buffer + offset;
	while (true) {
		value = strtol(start, &end, 10);
		if (start != end) {
			all_longs[i] = value;
			i++;
		} else
			break;
		start = end;
		end = NULL;
	}

	return i;
}
