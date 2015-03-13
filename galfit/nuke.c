#include <curses.h>

void printnuke ()
{
    extern char device[];
    extern int (*pfunc) (const char *, ...);

    pfunc ("                 __/~*##$%%@@@******~\\-__        \n");
    pfunc ("               /f=r/~_-~ _-_ --_.^-~--\\=b\\      \n");
    pfunc ("             4fF / */  .o  ._-__.__/~-. \\*R\\    \n");
    pfunc ("            /fF./  . /- /' /|/|  \\_  * *\\ *\\R\\  \n");
    pfunc ("           (iC.I+ '| - *-/00  |-  \\  )  ) )|RB  \n");
    pfunc ("           (I| (  [  / -|/^^\\ |   )  /_/ | *)B  \n");
    pfunc ("           (I(. \\ `` \\   \\m_m_|~__/ )_ .-~ F/   \n");
    pfunc ("            \\b\\\\=_.\\_b`-+-~x-_/ .. ,._/ , F/    \n");
    pfunc ("             ~\\_\\= =  =-*###%%#x==-#  *=- =/     \n");
    pfunc ("                ~\\**U/~  | i i | ~~~\\===~       \n");
    pfunc ("                        | I I \\\\                \n");
    pfunc ("                       / // i\\ \\\\               \n");
    pfunc ("                  (   [ (( I@) )))  )           \n");
    pfunc ("                       \\_\\_VYVU_/               \n");
    pfunc ("                         || * |                 \n");
    pfunc ("                        /* /I\\ *~~\\              \n");
    pfunc ("                      /~-/*  / \\ \\ ~~M~\\         \n");
    pfunc ("            ____----=~ // /WVW\\* \\|\\ ***===--___ \n");
    pfunc (" \n");
    if (strncmp (device, "regular", 7) != 0)
        refresh();
}


