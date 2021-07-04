#include <stdio.h>
#include <stdlib.h>

int main(){
    FILE *F1,*F2;
    F1=fopen("Messages out.txt", "r");
    if(F1==NULL) printf("Messages out - NULL");

    F2=fopen("Messages in.txt", "r");
    if(F2==NULL) printf("Messages in - NULL");

    int m_in, m_out,i=0;

    while((fscanf(F1,"%d", &(m_out))!=EOF)&&(fscanf(F2,"%d", &(m_in))!=EOF)){
        if(m_out!=m_in) printf("Error! %d\n", i);
        if(!(i%1000000)) printf("Symbols %d\n", i);
        i++;
    }
    return 0;
}
