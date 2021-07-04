#include <stdio.h>
#include <stdlib.h>

int Ht[24][12]={
    {1,0,0,0,0,0,0,0,0,0,0,0},
    {0,1,0,0,0,0,0,0,0,0,0,0},
    {0,0,1,0,0,0,0,0,0,0,0,0},
    {0,0,0,1,0,0,0,0,0,0,0,0},
    {0,0,0,0,1,0,0,0,0,0,0,0},
    {0,0,0,0,0,1,0,0,0,0,0,0},
    {0,0,0,0,0,0,1,0,0,0,0,0},
    {0,0,0,0,0,0,0,1,0,0,0,0},
    {0,0,0,0,0,0,0,0,1,0,0,0},
    {0,0,0,0,0,0,0,0,0,1,0,0},
    {0,0,0,0,0,0,0,0,0,0,1,0},
    {0,0,0,0,0,0,0,0,0,0,0,1},
    {1,1,0,1,1,1,0,0,0,1,0,1},
    {1,0,1,1,1,0,0,0,1,0,1,1},
    {0,1,1,1,0,0,0,1,0,1,1,1},
    {1,1,1,0,0,0,1,0,1,1,0,1},
    {1,1,0,0,0,1,0,1,1,0,1,1},
    {1,0,0,0,1,0,1,1,0,1,1,1},
    {0,0,0,1,0,1,1,0,1,1,1,1},
    {0,0,1,0,1,1,0,1,1,1,0,1},
    {0,1,0,1,1,0,1,1,1,0,0,1},
    {1,0,1,1,0,1,1,1,0,0,0,1},
    {0,1,1,0,1,1,1,0,0,0,1,1},
    {1,1,1,1,1,1,1,1,1,1,1,0}
};

int G[12][24]={
    {1,1,0,1,1,1,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0},
    {1,0,1,1,1,0,0,0,1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0},
    {0,1,1,1,0,0,0,1,0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0},
    {1,1,1,0,0,0,1,0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0},
    {1,1,0,0,0,1,0,1,1,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0},
    {1,0,0,0,1,0,1,1,0,1,1,1,0,0,0,0,0,1,0,0,0,0,0,0},
    {0,0,0,1,0,1,1,0,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0},
    {0,0,1,0,1,1,0,1,1,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0},
    {0,1,0,1,1,0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0},
    {1,0,1,1,0,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0},
    {0,1,1,0,1,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0},
    {1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1}
};

void zero (int mass[], int n){
    for(int i=0;i<n;i++) mass[i]=0;
}

void sum (int result[],int i){
    for(int j=0;j<12;j++){
        if((result[j])^(G[i][j])) result[j]=1;
        else result[j]=0;
    }
}


void Encoding (int u[], int result[]){   //"coder" gets mail (mass) of 12 bits and converts they in code word (mass) of 24 bits
    zero(result,24);
    for(int i=0;i<12;i++){
        result[i+12]=u[i];
        if(u[i]) sum(result,i);
    }
}

int main (void){

    int message[12]={0,0,0,0,0,0,0,0,0,0,0,0}, vector[24];

    FILE *F1,*F2;
    F1=fopen ("Messages out.txt", "w");
    if (F1==NULL) printf ("Messages out - NULL"), exit(2);

    F2=fopen ("Encoding with errors.txt", "w");
    if (F2==NULL) printf ("Encoding with errors - NULL"), exit(2);

    int a=0,b,c,y1,y2,y3;
    for(int i=0;i<4095;i++){
        for(int j=11;j>=0;j--){
            if(j==11) c=1;
            else c=0;
            b=a+message[j]+c;
            switch(b){
                case 3:message[j]=1, a=1;break;
                case 2:message[j]=0, a=1;break;
                case 1:message[j]=1, a=0;break;
                case 0:message[j]=0, a=0;break;
            }
        }
        a=0;

        Encoding(message,vector);

        for(int x1=0;x1<24;x1++){
            for(int x2=(x1+1);x2<23;x2++){
                for(int x3=(x2+1);x3<24;x3++){
                    y1=vector[x1];
                    y2=vector[x2];
                    y3=vector[x3];
                    if(vector[x1]) vector[x1]=0;
                    else vector[x1]=1;
                    if(vector[x2]) vector[x2]=0;
                    else vector[x2]=1;
                    if(vector[x3]) vector[x3]=0;
                    else vector[x3]=1;
                    for (int l=0;l<24;l++) fprintf (F2, "%d ",vector[l]);
                    fprintf (F2, "\n");
                    vector[x1]=y1;
                    vector[x2]=y2;
                    vector[x3]=y3;
                    for(int j=0;j<12;j++) fprintf(F1,"%d ", message[j]);
                    fprintf (F1, "\n");
                }
            }
        }
    }
    fclose(F1);
    fclose(F2);

    return 0;
}
