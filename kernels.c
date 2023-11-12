/********************************************************
 * Kernels to be optimized for the CS:APP Performance Lab
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"

/*
 * Please fill in the following team struct
 */
team_t team = {
    "pmat", /* Team name */

    "Patrick Bohn Matthiesen", /* First member full name */
    "pmat@itu.dk",             /* First member email address */

    "", /* Second member full name (leave blank if none) */
    ""  /* Second member email addr (leave blank if none) */
};

/***************
 * ROTATE KERNEL
 ***************/

/******************************************************
 * Your different versions of the rotate kernel go here
 ******************************************************/

/*
 * naive_rotate - The naive baseline version of rotate
 */
char naive_rotate_descr[] = "naive_rotate: Naive baseline implementation";
void naive_rotate(int dim, pixel *src, pixel *dst)
{
    int i, j;

    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
            dst[RIDX(dim - 1 - j, i, dim)] = src[RIDX(i, j, dim)];
}

/*
 * rotate - Your current working version of rotate
 * IMPORTANT: This is the version you will be graded on
 */
char rotate_descr[] = "rotate: Current working version, using loop unrolling and locality optimization";
void rotate(int dim, pixel *src, pixel *dst)
{
    int i, j;
    int dimMin1 = dim - 1;

    for (i = 0; i < dim; i += 16)
    {
        for (j = dimMin1; j >= 0; j--)
        {
            int jDim = RIDX(dimMin1 - j, i, dim);
            int iDim = RIDX(i, j, dim);

            dst[jDim++] = src[iDim];
            dst[jDim++] = src[iDim + dim];
            dst[jDim++] = src[iDim + 2 * dim];
            dst[jDim++] = src[iDim + 3 * dim];
            dst[jDim++] = src[iDim + 4 * dim];
            dst[jDim++] = src[iDim + 5 * dim];
            dst[jDim++] = src[iDim + 6 * dim];
            dst[jDim++] = src[iDim + 7 * dim];
            dst[jDim++] = src[iDim + 8 * dim];
            dst[jDim++] = src[iDim + 9 * dim];
            dst[jDim++] = src[iDim + 10 * dim];
            dst[jDim++] = src[iDim + 11 * dim];
            dst[jDim++] = src[iDim + 12 * dim];
            dst[jDim++] = src[iDim + 13 * dim];
            dst[jDim++] = src[iDim + 14 * dim];
            dst[jDim] = src[iDim + 15 * dim];
        }
    }
}

char rotate_descr4[] = "rotate: Loop has been unrolled";
void rotate4(int dim, pixel *src, pixel *dst)
{
    int i, j;
    for (i = 0; i < dim; i += 16)
        for (j = 0; j < dim; j++)
        {
            dst[RIDX(dim - 1 - j, i + 0, dim)] = src[RIDX(i + 0, j, dim)];
            dst[RIDX(dim - 1 - j, i + 1, dim)] = src[RIDX(i + 1, j, dim)];
            dst[RIDX(dim - 1 - j, i + 2, dim)] = src[RIDX(i + 2, j, dim)];
            dst[RIDX(dim - 1 - j, i + 3, dim)] = src[RIDX(i + 3, j, dim)];
            dst[RIDX(dim - 1 - j, i + 4, dim)] = src[RIDX(i + 4, j, dim)];
            dst[RIDX(dim - 1 - j, i + 5, dim)] = src[RIDX(i + 5, j, dim)];
            dst[RIDX(dim - 1 - j, i + 6, dim)] = src[RIDX(i + 6, j, dim)];
            dst[RIDX(dim - 1 - j, i + 7, dim)] = src[RIDX(i + 7, j, dim)];
            dst[RIDX(dim - 1 - j, i + 8, dim)] = src[RIDX(i + 8, j, dim)];
            dst[RIDX(dim - 1 - j, i + 9, dim)] = src[RIDX(i + 9, j, dim)];
            dst[RIDX(dim - 1 - j, i + 10, dim)] = src[RIDX(i + 10, j, dim)];
            dst[RIDX(dim - 1 - j, i + 11, dim)] = src[RIDX(i + 11, j, dim)];
            dst[RIDX(dim - 1 - j, i + 12, dim)] = src[RIDX(i + 12, j, dim)];
            dst[RIDX(dim - 1 - j, i + 13, dim)] = src[RIDX(i + 13, j, dim)];
            dst[RIDX(dim - 1 - j, i + 14, dim)] = src[RIDX(i + 14, j, dim)];
            dst[RIDX(dim - 1 - j, i + 15, dim)] = src[RIDX(i + 15, j, dim)];
        }
}

/*********************************************************************
 * register_rotate_functions - Register all of your different versions
 *     of the rotate kernel with the driver by calling the
 *     add_rotate_function() for each test function. When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.
 *********************************************************************/

void register_rotate_functions()
{
    add_rotate_function(&naive_rotate, naive_rotate_descr);
    add_rotate_function(&rotate, rotate_descr);
    add_rotate_function(&rotate4, rotate_descr4);
    /* ... Register additional test functions here */
}

/***************
 * SMOOTH KERNEL
 **************/

/***************************************************************
 * Various typedefs and helper functions for the smooth function
 * You may modify these any way you like.
 **************************************************************/

/* A struct used to compute averaged pixel value */
typedef struct
{
    int red;
    int green;
    int blue;
    int num;
} pixel_sum;

/* Compute min and max of two integers, respectively */
static int min(int a, int b) { return (a < b ? a : b); }
static int max(int a, int b) { return (a > b ? a : b); }

/*
 * initialize_pixel_sum - Initializes all fields of sum to 0
 */
static void initialize_pixel_sum(pixel_sum *sum)
{
    sum->red = sum->green = sum->blue = 0;
    sum->num = 0;
    return;
}

/*
 * accumulate_sum - Accumulates field values of p in corresponding
 * fields of sum
 */
static void accumulate_sum(pixel_sum *sum, pixel p)
{
    sum->red += (int)p.red;
    sum->green += (int)p.green;
    sum->blue += (int)p.blue;
    sum->num++;
    return;
}

/*
 * assign_sum_to_pixel - Computes averaged pixel value in current_pixel
 */
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum)
{
    current_pixel->red = (unsigned short)(sum.red / sum.num);
    current_pixel->green = (unsigned short)(sum.green / sum.num);
    current_pixel->blue = (unsigned short)(sum.blue / sum.num);
    return;
}

/*
 * avg - Returns averaged pixel value at (i,j)
 */
static pixel avg(int dim, int i, int j, pixel *src)
{
    int ii, jj;
    pixel_sum sum; // extract to 3 sums, one for each line, it should give better use of the pipelines
    pixel current_pixel;

    initialize_pixel_sum(&sum);

    int iMin = min(i + 1, dim - 1);
    int jMin = min(j + 1, dim - 1);

    for (ii = max(i - 1, 0); ii <= iMin; ii++)
        for (jj = max(j - 1, 0); jj <= jMin; jj++)
            accumulate_sum(&sum, src[RIDX(ii, jj, dim)]);

    assign_sum_to_pixel(&current_pixel, sum);
    return current_pixel;
}

/******************************************************
 * Your different versions of the smooth kernel go here
 ******************************************************/

/*
 * naive_smooth - The naive baseline version of smooth
 */
char naive_smooth_descr[] = "naive_smooth: Naive baseline implementation";
void naive_smooth(int dim, pixel *src, pixel *dst)
{
    int i, j;

    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
            dst[RIDX(i, j, dim)] = avg(dim, i, j, src);
}

static inline void cornor(int index, int inc, int dim, pixel *src, pixel *dst)
{
    dst[index].red = (src[index].red + src[index + inc].red + src[index + dim].red + src[index + inc + dim].red) >> 2;
    dst[index].green = (src[index].green + src[index + inc].green + src[index + dim].green + src[index + inc + dim].green) >> 2;
    dst[index].blue = (src[index].blue + src[index + inc].blue + src[index + dim].blue + src[index + inc + dim].blue) >> 2;
}

static inline void rowEdge(int index, int dim, pixel *src, unsigned int *r, unsigned int *g, unsigned int *b)
{
    *r = src[index].red + src[index + 1].red + src[index + 2].red + src[index + dim].red + src[index + 1 + dim].red + src[index + 2 + dim].red;
    *g = src[index].green + src[index + 1].green + src[index + 2].green + src[index + dim].green + src[index + 1 + dim].green + src[index + 2 + dim].green;
    *b = src[index].blue + src[index + 1].blue + src[index + 2].blue + src[index + dim].blue + src[index + 1 + dim].blue + src[index + 2 + dim].blue;
}

static inline void colEdge(int index, int dim, pixel *src, unsigned int *r, unsigned int *g, unsigned int *b)
{
    *r = src[index].red + src[index + dim].red + src[index + (dim << 1)].red + src[index + 1].red + src[index + 1 + dim].red + src[index + 1 + (dim << 1)].red;
    *g = src[index].green + src[index + dim].green + src[index + (dim << 1)].green + src[index + 1].green + src[index + 1 + dim].green + src[index + 1 + (dim << 1)].green;
    *b = src[index].blue + src[index + dim].blue + src[index + (dim << 1)].blue + src[index + 1].blue + src[index + 1 + dim].blue + src[index + 1 + (dim << 1)].blue;
}

static inline void moveToNext(int index, int inc, int index2, int inc2, pixel *src, unsigned int *r, unsigned int *g, unsigned int *b)
{
    *r += src[index].red + src[index + inc].red - src[index2].red - src[index2 + inc2].red;
    *g += src[index].green + src[index + inc].green - src[index2].green - src[index2 + inc2].green;
    *b += src[index].blue + src[index + inc].blue - src[index2].blue - src[index2 + inc2].blue;
}

static inline void addRowsX3(int index, int dim, pixel *src, unsigned int *r, unsigned int *g, unsigned int *b)
{
    *r += src[index].red + src[index + dim].red + src[index + dim + dim].red;
    *g += src[index].green + src[index + dim].green + src[index + dim + dim].green;
    *b += src[index].blue + src[index + dim].blue + src[index + dim + dim].blue;
}

static inline void removeRowsX3(int index, int dim, pixel *src, unsigned int *r, unsigned int *g, unsigned int *b)
{
    *r -= src[index].red + src[index + dim].red + src[index + dim + dim].red;
    *g -= src[index].green + src[index + dim].green + src[index + dim + dim].green;
    *b -= src[index].blue + src[index + dim].blue + src[index + dim + dim].blue;
}

char smooth_descr[] = "smooth: Current working version, using blocking and saving of old pixel sums";
void smooth(int dim, pixel *src, pixel *dst)
{
    int i, j, k;

    int blocksize = 32; // 16 is slightly faster on average, but it has worse performance as the image size increases
    unsigned int sums[(blocksize + 2) * 3 + 3];
    unsigned int sumr, sumg, sumb; // i hope that it adds this to the register as it is used a lot

    // for each block then sum 3 pixels under each other in a collumn, and do that for the block length
    // adding together the sum of 3 columns of 3 pixels (9 in total), setting the middle pixel to the average of the 9 pixels
    // then move the block down 1 pixel and repeat.
    // it runs on blocksize amount of columns at a time, from row 1 to row dim - 1
    for (j = 0; j < dim; j += blocksize)
    {
        // get initial sums for the first row
        for (k = 0; k < blocksize + 2; k++)
        {
            int index = RIDX(0, j + k, dim);
            sums[k * 3 + 0] = src[index].red + src[index + dim].red + src[index + (dim << 1)].red;
            sums[k * 3 + 1] = src[index].green + src[index + dim].green + src[index + (dim << 1)].green;
            sums[k * 3 + 2] = src[index].blue + src[index + dim].blue + src[index + (dim << 1)].blue;
        }
        // run from row 1 to row dim-2
        for (i = 1; i < dim - 1; i++)
        {
            // adding together the sum of 3 columns of 3 pixels (9 in total)
            sumr = sums[0] + sums[3] + sums[6];
            sumg = sums[1] + sums[4] + sums[7];
            sumb = sums[2] + sums[5] + sums[8];
            // for each pixel on the current row
            // set the middle pixel to the average of the 9 pixels
            // and update the sums by removeing the column that is to the left and add the column that is to the right
            for (k = j + 1; k < min(j + 1 + blocksize, dim); k++)
            {
                dst[RIDX(i, k, dim)].red = sumr / 9;
                dst[RIDX(i, k, dim)].green = sumg / 9;
                dst[RIDX(i, k, dim)].blue = sumb / 9;
                sumr += sums[(k - j - 1 + 3) * 3 + 0] - sums[(k - j - 1) * 3 + 0];
                sumg += sums[(k - j - 1 + 3) * 3 + 1] - sums[(k - j - 1) * 3 + 1];
                sumb += sums[(k - j - 1 + 3) * 3 + 2] - sums[(k - j - 1) * 3 + 2];
            }
            // move one line down, remove line above, and add new line under the current one
            for (k = 0; k < blocksize + 2; k++)
            {
                sums[k * 3 + 0] += src[RIDX(i + 2, j + k, dim)].red - src[RIDX(i - 1, j + k, dim)].red;
                sums[k * 3 + 1] += src[RIDX(i + 2, j + k, dim)].green - src[RIDX(i - 1, j + k, dim)].green;
                sums[k * 3 + 2] += src[RIDX(i + 2, j + k, dim)].blue - src[RIDX(i - 1, j + k, dim)].blue;
            }
        }
    }

    // for (j = 0; j < dim; j += blocksize)
    // {
    //     // run from row 1 to row dim-2
    //     for (i = 1; i < dim - 1; i++)
    //     {
    //         // adding together the sum of 3 columns of 3 pixels (9 in total)
    //         sumr = src[RIDX(i - 1, j, dim)].red + src[RIDX(i - 1, j + 1, dim)].red + src[RIDX(i - 1, j + 2, dim)].red + src[RIDX(i, j, dim)].red + src[RIDX(i, j + 1, dim)].red + src[RIDX(i, j + 2, dim)].red + src[RIDX(i + 1, j, dim)].red + src[RIDX(i + 1, j + 1, dim)].red + src[RIDX(i + 1, j + 2, dim)].red;
    //         sumg = src[RIDX(i - 1, j, dim)].green + src[RIDX(i - 1, j + 1, dim)].green + src[RIDX(i - 1, j + 2, dim)].green + src[RIDX(i, j, dim)].green + src[RIDX(i, j + 1, dim)].green + src[RIDX(i, j + 2, dim)].green + src[RIDX(i + 1, j, dim)].green + src[RIDX(i + 1, j + 1, dim)].green + src[RIDX(i + 1, j + 2, dim)].green;
    //         sumb = src[RIDX(i - 1, j, dim)].blue + src[RIDX(i - 1, j + 1, dim)].blue + src[RIDX(i - 1, j + 2, dim)].blue + src[RIDX(i, j, dim)].blue + src[RIDX(i, j + 1, dim)].blue + src[RIDX(i, j + 2, dim)].blue + src[RIDX(i + 1, j, dim)].blue + src[RIDX(i + 1, j + 1, dim)].blue + src[RIDX(i + 1, j + 2, dim)].blue;
    //         // for each pixel on the current row
    //         // set the middle pixel to the average of the 9 pixels
    //         // and update the sum by removeing the column that is to the left and add the column that is to the right
    //         for (k = j + 1; k < min(j + blocksize + 1, dim); k++) {
    //             dst[RIDX(i, k, dim)].red = sumr / 9;
    //             dst[RIDX(i, k, dim)].green = sumg / 9;
    //             dst[RIDX(i, k, dim)].blue = sumb / 9;
    //             sumr += src[RIDX(i - 1, k + 2, dim)].red + src[RIDX(i, k + 2, dim)].red + src[RIDX(i + 1, k + 2, dim)].red - src[RIDX(i - 1, k - 1, dim)].red - src[RIDX(i, k - 1, dim)].red - src[RIDX(i + 1, k - 1, dim)].red;
    //             sumg += src[RIDX(i - 1, k + 2, dim)].green + src[RIDX(i, k + 2, dim)].green + src[RIDX(i + 1, k + 2, dim)].green - src[RIDX(i - 1, k - 1, dim)].green - src[RIDX(i, k - 1, dim)].green - src[RIDX(i + 1, k - 1, dim)].green;
    //             sumb += src[RIDX(i - 1, k + 2, dim)].blue + src[RIDX(i, k + 2, dim)].blue + src[RIDX(i + 1, k + 2, dim)].blue - src[RIDX(i - 1, k - 1, dim)].blue - src[RIDX(i, k - 1, dim)].blue - src[RIDX(i + 1, k - 1, dim)].blue;
    //         }
    //     }
    // }

    // for each edge then take the initial sum of 6 pixels
    // then set the average of the 6 pixels in dst
    // then update the sum to reflect the new pixels that are going into and out of proximity of the next pixel

    // first row
    rowEdge(RIDX(0, 0, dim), dim, src, &sumr, &sumg, &sumb);
    for (i = 1; i < dim - 1; i++)
    {
        dst[RIDX(0, i, dim)].red = sumr / 6;
        dst[RIDX(0, i, dim)].green = sumg / 6;
        dst[RIDX(0, i, dim)].blue = sumb / 6;
        moveToNext(RIDX(0, i + 2, dim), dim, RIDX(0, i - 1, dim), dim, src, &sumr, &sumg, &sumb);
    }
    // last row
    rowEdge(RIDX(dim - 2, 0, dim), dim, src, &sumr, &sumg, &sumb);
    for (i = 1; i < dim - 1; i++)
    {
        dst[RIDX(dim - 1, i, dim)].red = sumr / 6;
        dst[RIDX(dim - 1, i, dim)].green = sumg / 6;
        dst[RIDX(dim - 1, i, dim)].blue = sumb / 6;
        moveToNext(RIDX(dim - 2, i + 2, dim), dim, RIDX(dim - 2, i - 1, dim), dim, src, &sumr, &sumg, &sumb);
    }

    // first col
    colEdge(RIDX(0, 0, dim), dim, src, &sumr, &sumg, &sumb);
    for (i = 1; i < dim - 1; i++)
    {
        dst[RIDX(i, 0, dim)].red = sumr / 6;
        dst[RIDX(i, 0, dim)].green = sumg / 6;
        dst[RIDX(i, 0, dim)].blue = sumb / 6;
        moveToNext(RIDX(i + 2, 0, dim), 1, RIDX(i - 1, 0, dim), 1, src, &sumr, &sumg, &sumb);
    }
    // last col
    colEdge(RIDX(0, dim - 2, dim), dim, src, &sumr, &sumg, &sumb);
    for (i = 1; i < dim - 1; i++)
    {
        dst[RIDX(i, dim - 1, dim)].red = sumr / 6;
        dst[RIDX(i, dim - 1, dim)].green = sumg / 6;
        dst[RIDX(i, dim - 1, dim)].blue = sumb / 6;
        moveToNext(RIDX(i + 2, dim - 2, dim), 1, RIDX(i - 1, dim - 2, dim), 1, src, &sumr, &sumg, &sumb);
    }

    // corners
    cornor(RIDX(0, 0, dim), 1, dim, src, dst);               // top left
    cornor(RIDX(0, dim - 1, dim), -1, dim, src, dst);        // top right
    cornor(RIDX(dim - 1, 0, dim), 1, -dim, src, dst);        // bottom left
    cornor(RIDX(dim - 1, dim - 1, dim), -1, -dim, src, dst); // bottom right
}

/*********************************************************************
 * register_smooth_functions - Register all of your different versions
 *     of the smooth kernel with the driver by calling the
 *     add_smooth_function() for each test function.  When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.
 *********************************************************************/

void register_smooth_functions()
{
    add_smooth_function(&naive_smooth, naive_smooth_descr);
    add_smooth_function(&smooth, smooth_descr);
    /* ... Register additional test functions here */
}
