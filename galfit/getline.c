/*  This subroutine reads key strokes into a variable called "buffer."
    Buffer will store at most n-1 characters, while the end of line
    character is tacked onto the nth character to terminate the string.
    Strings longer than the allowed number of spaces are read, but
    discarded.  This is to prevent extra characters from being automatically
    interpreted by the next "scanf" call, which can often produce undesired
    results.								*/

# include <stdio.h>

void getline (char buffer[], int n)
{
    char character;
    int i=0;

    do {
	character = getchar();
	if (i < n) buffer[i++] = character;
    } while (character != '\n');
    buffer[i-1] = '\0';
}
