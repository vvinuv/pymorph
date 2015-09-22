
/* This subroutine clears the keyboard buffer.   This can be used after
   a "scanf" statement to get rid of extra characters that the user 
   accidentally entered, which otherwise may screw up the next "scanf" 
   command 								*/

void getline(char *, int);

void clearbuff(void) 
{
    char cbuff[2];

    getline(cbuff, 2);
}
