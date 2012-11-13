#include <stdio.h>
#include <stdbool.h>

bool isFloat(char* c){
  while (1){
    if((*c)=='.')
      return true;
    else if((*c)=='\0')
      return false;
    c++;
  }
}

int main(int argc, char* argv[]){
  if( argc < 3 ){
    printf("Usage:\nbinconv <input file> <output file>\n");
    return -1;
  }
  FILE *input_file = fopen(argv[1],"r");
  FILE *output_file = fopen(argv[2],"wb"); 
  double g;
  int i;
  char line[50];
  while( fscanf(input_file,"%s", line ) != EOF  ){    
    if( !isFloat(line )){
      sscanf(line,"%d",&i);
      fwrite(&i,sizeof(int),1,output_file);
    }
    else{
      sscanf(line,"%lg",&g);
      fwrite(&g,sizeof(double),1,output_file);
    }
  }
  fclose(input_file);
  fclose(output_file); 
  return 0;
}
