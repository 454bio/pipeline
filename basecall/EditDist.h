#pragma once

// distance - performs the basic edit distance with penalties to convert the second word into the first
// note: the first word is the reference, and the second word is the string to edit in order to match the reference
int distance(const char *word1, int len1, const char *word2, int len2);

