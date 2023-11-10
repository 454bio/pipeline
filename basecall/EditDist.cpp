#include <stdio.h>
#include <string.h>

int penalty_delete = 1;
int penalty_insert = 1;
int penalty_substitute = 1;

// distance - performs the basic edit distance with penalties to convert the second word into the first
// note: the first word is the reference, and the second word is the string to edit in order to match the reference
int distance(const char *word1, int len1, const char *word2, int len2)
{
    int matrix[len1 + 1][len2 + 1];
    int i;
    for (i = 0; i <= len1; i++) {
        matrix[i][0] = i;
    }
    for (i = 0; i <= len2; i++) {
        matrix[0][i] = i;
    }
    for (i = 1; i <= len1; i++) {
        int j;
        char c1;

        c1 = word1[i-1];
        for (j = 1; j <= len2; j++) {
            char c2;

            c2 = word2[j-1];
            if (c1 == c2) {
                matrix[i][j] = matrix[i-1][j-1];
            } else {
                int cost_delete;
                int cost_insert;
                int cost_substitute;
                int cost_minimum;

                cost_delete = matrix[i-1][j] + penalty_delete;
                cost_insert = matrix[i][j-1] + penalty_insert;
                cost_substitute = matrix[i-1][j-1] + penalty_substitute;
                cost_minimum = cost_delete;
                if (cost_insert < cost_minimum) {
                    cost_minimum = cost_insert;
                }
                if (cost_substitute < cost_minimum) {
                    cost_minimum = cost_substitute;
                }
                matrix[i][j] = cost_minimum;
            }
        }
    }
    return matrix[len1][len2];
}

