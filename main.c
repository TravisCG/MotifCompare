#include <stdlib.h>
#include <stdio.h>

#include "motif.h"

int main(int argc, char **argv){
   Motif motif;

   if(argc < 2){
      printf("Please give me the search directory (do not forget the slash!)\n");
      printf("Example: ./moco /home/travis/fancymotifs/\n");
      return(1);
   }

   motif = readDir(argv[1]);
   motifAllAgainstAll(motif);

   return(EXIT_SUCCESS);
}
