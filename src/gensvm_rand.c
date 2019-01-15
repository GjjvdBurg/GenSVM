/**
 * @file gensvm_rand.c
 * @author G.J.J. van den Burg
 * @date 2018-04-17
 * @brief Functions for random number generation
 *
 * @copyright
 * Copyright 2018, G.J.J. van den Burg

 GenSVM is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 GenSVM is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with GenSVM; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

 */

#include "gensvm_rand.h"

/* Linear congruential.  */
#define	TYPE_0		0
#define	BREAK_0		8
#define	DEG_0		0
#define	SEP_0		0

/* x**7 + x**3 + 1.  */
#define	TYPE_1		1
#define	BREAK_1		32
#define	DEG_1		7
#define	SEP_1		3

/* x**15 + x + 1.  */
#define	TYPE_2		2
#define	BREAK_2		64
#define	DEG_2		15
#define	SEP_2		1

/* x**31 + x**3 + 1.  */
#define	TYPE_3		3
#define	BREAK_3		128
#define	DEG_3		31
#define	SEP_3		3

/* x**63 + x + 1.  */
#define	TYPE_4		4
#define	BREAK_4		256
#define	DEG_4		63
#define	SEP_4		1


/* Array versions of the above information to make code run faster.
   Relies on fact that TYPE_i == i.  */

#define	MAX_TYPES	5	/* Max number of types above.  */

struct random_poly_info
{
	int seps[MAX_TYPES];
	int degrees[MAX_TYPES];
};

static const struct random_poly_info random_poly_info = 
{
	{ SEP_0, SEP_1, SEP_2, SEP_3, SEP_4 },
	{ DEG_0, DEG_1, DEG_2, DEG_3, DEG_4 }
};

static int32_t randtbl[DEG_3 + 1] =
{
	TYPE_3,

	-1726662223, 379960547, 1735697613, 1040273694, 1313901226,
	1627687941, -179304937, -2073333483, 1780058412, -1989503057,
	-812545463, 154635395, 1388815473, -1926676823, 525320961,
	-1009028674, 968117788, -123449607, 1284210865, 435012392,
	-2017506339, -911064859, -370259173, 1132637927, 1398500161,
	-205601318,
};

static struct gensvm_random_data unsafe_state =
{
	.fptr = &randtbl[SEP_3 + 1],
	.rptr = &randtbl[1],
	.state = &randtbl[1],
	.rand_type = TYPE_3,
	.rand_deg = DEG_3,
	.rand_sep = SEP_3,
	.end_ptr = &randtbl[sizeof (randtbl) / sizeof(randtbl[0])]
};

char *gensvm_rand_initstate(unsigned int seed, char *arg_state, size_t n);
char *gensvm_rand_setstate(char *arg_state);
int gensvm_rand_srandom_r (unsigned int seed, struct gensvm_random_data *buf);
int gensvm_rand_initstate_r (unsigned int seed, char *arg_state, size_t n, struct gensvm_random_data *buf);
int gensvm_rand_setstate_r (char *arg_state, struct gensvm_random_data *buf);
int gensvm_rand_random_r (struct gensvm_random_data *buf, int32_t *result);
long int gensvm_random(void);
void gensvm_srandom(unsigned int x);

/** @brief Replacement for builtin rand()
 *
 * @return random number
 */
int gensvm_rand()
{
	return (int) gensvm_random();
}

/** @brief Replacement for builtin srand()
 *
 * @param x seed
 */
void gensvm_srand(unsigned int x)
{
	(void) gensvm_srandom(x);
}

void gensvm_srandom(unsigned int x)
{
	(void) gensvm_rand_srandom_r(x, &unsafe_state);
}

char *gensvm_rand_initstate(unsigned int seed, char *arg_state, size_t n)
{
	int32_t *ostate;

	ostate = &unsafe_state.state[-1];
	gensvm_rand_initstate_r(seed, arg_state, n, &unsafe_state);

	return (char *) ostate;
}

char *gensvm_rand_setstate(char *arg_state)
{
	int32_t *ostate;

	ostate = &unsafe_state.state[-1];

	if (gensvm_rand_setstate_r (arg_state, &unsafe_state) < 0)
		ostate = NULL;

	return (char *) ostate;
}

long int gensvm_random()
{
	int32_t retval;

	(void) gensvm_rand_random_r(&unsafe_state, &retval);

	return retval;
}

int gensvm_rand_srandom_r (unsigned int seed, struct gensvm_random_data *buf)
{
	int type;
	int32_t *state;
	long int i;
	int32_t word;
	int32_t *dst;
	int kc;

	if (buf == NULL)
		goto fail;
	type = buf->rand_type;
	if ((unsigned int) type >= MAX_TYPES)
		goto fail;

	state = buf->state;
	/* We must make sure the seed is not 0.  Take arbitrarily 1 in this case.  */
	if (seed == 0)
		seed = 1;
	state[0] = seed;
	if (type == TYPE_0)
		goto done;

	dst = state;
	word = seed;
	kc = buf->rand_deg;
	for (i = 1; i < kc; ++i)
	{
		/* This does:
		   state[i] = (16807 * state[i - 1]) % 2147483647;
		   but avoids overflowing 31 bits.  */
		long int hi = word / 127773;
		long int lo = word % 127773;
		word = 16807 * lo - 2836 * hi;
		if (word < 0)
			word += 2147483647;
		*++dst = word;
	}

	buf->fptr = &state[buf->rand_sep];
	buf->rptr = &state[0];
	kc *= 10;
	while (--kc >= 0)
	{
		int32_t discard;
		(void) gensvm_rand_random_r (buf, &discard);
	}

done:
	return 0;

fail:
	return -1;
}

int gensvm_rand_initstate_r (unsigned int seed, char *arg_state, size_t n, struct gensvm_random_data *buf)
{
	if (buf == NULL)
		goto fail;

	int32_t *old_state = buf->state;
	if (old_state != NULL)
	{
		int old_type = buf->rand_type;
		if (old_type == TYPE_0)
			old_state[-1] = TYPE_0;
		else
			old_state[-1] = (MAX_TYPES * (buf->rptr - old_state)) + old_type;
	}

	int type;
	if (n >= BREAK_3)
		type = n < BREAK_4 ? TYPE_3 : TYPE_4;
	else if (n < BREAK_1)
	{
		if (n < BREAK_0)
			goto fail;

		type = TYPE_0;
	}
	else
		type = n < BREAK_2 ? TYPE_1 : TYPE_2;

	int degree = random_poly_info.degrees[type];
	int separation = random_poly_info.seps[type];

	buf->rand_type = type;
	buf->rand_sep = separation;
	buf->rand_deg = degree;
	int32_t *state = &((int32_t *) arg_state)[1];	/* First location.  */
	/* Must set END_PTR before srandom.  */
	buf->end_ptr = &state[degree];

	buf->state = state;

	gensvm_rand_srandom_r (seed, buf);

	state[-1] = TYPE_0;
	if (type != TYPE_0)
		state[-1] = (buf->rptr - state) * MAX_TYPES + type;

	return 0;

fail:
	errno = EINVAL;
	return -1;
}


int gensvm_rand_setstate_r (char *arg_state, struct gensvm_random_data *buf)
{
	int32_t *new_state = 1 + (int32_t *) arg_state;
	int type;
	int old_type;
	int32_t *old_state;
	int degree;
	int separation;

	if (arg_state == NULL || buf == NULL)
		goto fail;

	old_type = buf->rand_type;
	old_state = buf->state;
	if (old_type == TYPE_0)
		old_state[-1] = TYPE_0;
	else
		old_state[-1] = (MAX_TYPES * (buf->rptr - old_state)) + old_type;

	type = new_state[-1] % MAX_TYPES;
	if (type < TYPE_0 || type > TYPE_4)
		goto fail;

	buf->rand_deg = degree = random_poly_info.degrees[type];
	buf->rand_sep = separation = random_poly_info.seps[type];
	buf->rand_type = type;

	if (type != TYPE_0)
	{
		int rear = new_state[-1] / MAX_TYPES;
		buf->rptr = &new_state[rear];
		buf->fptr = &new_state[(rear + separation) % degree];
	}
	buf->state = new_state;
	/* Set end_ptr too.  */
	buf->end_ptr = &new_state[degree];

	return 0;

fail:
	errno = EINVAL;
	return -1;
}

// Stdlib's rand() relies on integer overflow in the expression:
//
//     val = *fptr += *rptr;
//
// Sometimes packages are checked using the gcc undefined behavior sanitizer, 
// which throws a runtime error when signed integer overflow occurs.
//
// Because we want to stay close to the stdlib implementation of rand(), we 
// use this function to compute the sum of two integers that can overflow, 
// without actually allowing it to overflow.
int32_t overflow_add(int32_t a, int32_t b)
{
	int32_t p, q, res;
	if (a > 0 && b > INT32_MAX - a) {
		p = a - INT32_MAX;
		p = overflow_add(p, b);
		res = INT32_MIN + p - 1;
	} else if (a < 0 && b < INT32_MIN - a) {
		q = a + INT32_MAX;
		q = overflow_add(q, b);
		q = overflow_add(q, INT32_MAX);
		res = overflow_add(q, 2);
	} else {
		res = a + b;
	}
	return res;
}

int gensvm_rand_random_r (struct gensvm_random_data *buf, int32_t *result)
{
	int32_t *state;

	if (buf == NULL || result == NULL)
		goto fail;

	state = buf->state;

	if (buf->rand_type == TYPE_0)
	{
		int32_t val = state[0];
		val = ((state[0] * 1103515245) + 12345) & 0x7fffffff;
		state[0] = val;
		*result = val;
	}
	else
	{
		int32_t *fptr = buf->fptr;
		int32_t *rptr = buf->rptr;
		int32_t *end_ptr = buf->end_ptr;
		int32_t val;

		//val = *fptr += *rptr;
		*fptr = overflow_add(*fptr, *rptr);
		val = *fptr;

		/* Chucking least random bit.  */
		*result = (val >> 1) & 0x7fffffff;
		++fptr;
		if (fptr >= end_ptr)
		{
			fptr = state;
			++rptr;
		}
		else
		{
			++rptr;
			if (rptr >= end_ptr)
				rptr = state;
		}
		buf->fptr = fptr;
		buf->rptr = rptr;
	}
	return 0;

fail:
	errno = EINVAL;
	return -1;
}
