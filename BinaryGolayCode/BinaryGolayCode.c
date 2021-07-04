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

void copy (int mass1[], int mass2[], int n, int s){
    for(int i=0;i<n;i++) mass1[i]=mass2[i+s];
}

void binary_sum(int mass[],int SwapVector[],int message[]){
    for(int i=0;i<12;i++) {
        if(mass[i]^SwapVector[i]) message[i]=1;
        else message[i]=0;
    }
}

void sum (int result[],int i){
    for(int j=0;j<12;j++){
        if((result[j])^(G[i][j])) result[j]=1;
        else result[j]=0;
    }
}

int weight (int Syndrom[]){
    int W=0;
    for(int i=0;i<12;i++) W+=Syndrom[i];
    return W;
}

int number (int Syndrom[]){
    int mass[12];
    copy(mass,Syndrom,12,0);
    for(int i=0;i<12;i++){
            sum(mass,i);
            if((weight(mass))<=2) return i;
            else copy(mass,Syndrom,12,0);
    }
    return -1;
}

void Encoding (int u[], int result[]){
    zero(result,24);
    for(int i=0;i<12;i++){
        result[i+12]=u[i];
        if(u[i]) sum(result,i);
    }
}

void Syndrome (int CodeVector[], int Syndrom[]){
    zero(Syndrom,12);
    for(int j=0;j<12;j++){
        for(int i=0;i<24;i++){
            if((CodeVector[i])&&(Ht[i][j])) Syndrom[j]+=1;
        }
        if(!(Syndrom[j]%2)) Syndrom[j]=0;
        else Syndrom[j] = 1;
    }
}

void Decoding (int CodeVector[], int message[]){
    int Syndrom[12];
    zero(message,12);
    Syndrome(CodeVector,Syndrom);
    if((weight(Syndrom))<=3) copy(message,CodeVector,12,12);
    else{
        int mass[12],n;
        copy(mass,Syndrom,12,0);
        n=number(Syndrom);
        if(n>=0){
            copy(message,CodeVector,12,12);
            if(message[n])message[n]=0;
            else message[n]=1;
        }
        else {
            int temp, SwapVector[24];
            copy(SwapVector,CodeVector,24,0);
            for(int i=0;i<12;i++){
                temp=SwapVector[i];
                SwapVector[i]=SwapVector[i+12];
                SwapVector[i+12]=temp;
            }
            Syndrome(SwapVector,Syndrom);
            if((weight(Syndrom))<=3) binary_sum(Syndrom,SwapVector,message);
            else{
                copy(mass,Syndrom,12,0);
                n=number(Syndrom);
                sum(mass,n);
                if(n>=0) binary_sum(mass,SwapVector,message);
                else printf("Request retransmission.\n");
            }
        }
    }
}
