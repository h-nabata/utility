#include <stdio.h>
#include <math.h>

// 1~1000 までの 5n±1 型の整数、平方根の総和を計算するプログラム

int main(void) {
    FILE *fp;
    int sum_integer = 0;
    double sum_float = 0.0;
    int i;
    char fname[10] = "sum.txt";

    fp = fopen(fname,"w"); // open the file for writing
    if (fp == NULL) { // when failed to open
        printf("failed to open %s!\n", fname);
        return -1;
    }

    i = 1; // initialize i
    while (i <= 1000) {
        if (i % 5 == 1 || i % 5 == 4){
            sum_integer = sum_integer + i;
            fprintf(fp, "%2d\n", i); // output to the file
        }
        i++;
    }

    for (i = 1; i <= 1000; i++) {
        if (i % 5 == 1 || i % 5 == 4){
            sum_float = sum_float + sqrt(i);
            fprintf(fp, "%2d :\t%f\n", i, sqrt(i)); // output to the file
        }
    }

    printf(" int sum = %12d\n", sum_integer);
    printf("sqrt sum = %12.5f\n", sum_float);

    fclose(fp); // close the file
    return 0;
}
