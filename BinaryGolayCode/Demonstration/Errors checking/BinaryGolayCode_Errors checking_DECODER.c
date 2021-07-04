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

    //=>we have 0-3 errors in parity-check part (or don't have errors in all vector)
    if((weight(Syndrom))<=3) copy(message,CodeVector,12,12);

    else{//assumption: we have only one error in massage and 0-2 errors in parity-check part
        int mass[12],n;
        copy(mass,Syndrom,12,0);
        n=number(Syndrom);
        if(n>=0){//=>we have only one error in message and 0-2 errors in parity-check part
            copy(message,CodeVector,12,12);
            if(message[n])message[n]=0;
            else message[n]=1;
        }
        else {//assumption: we have 2-3 errors in message and 0-1 errors in parity-check part
            int temp, SwapVector[24];
            copy(SwapVector,CodeVector,24,0);
            for(int i=0;i<12;i++){//shifted left by 12 bits (swap message and parity-check part in vector)
                temp=SwapVector[i];
                SwapVector[i]=SwapVector[i+12];
                SwapVector[i+12]=temp;
            }
            Syndrome(SwapVector,Syndrom);

            //=> we have 2-3 errors only in message (0 error in parity-check part)
            if((weight(Syndrom))<=3) binary_sum(Syndrom,SwapVector,message);
            else{//assumption: we have 2 errors in message and 1 error in parity-check part
                copy(mass,Syndrom,12,0);
                n=number(Syndrom);
                sum(mass,n);

                //we have 2 errors in message and 1 error in parity-check part
                if(n>=0) binary_sum(mass,SwapVector,message);
                else printf("Request retransmission.\n");//we have 4 or more errors.
            }
        }
    }
}



int main (void){
    FILE *F1,*F2;

    F1=fopen ("Encoding with errors.txt", "r");
    if (F1==NULL) printf ("Encoding with errors - NULL"), exit(2);

    F2=fopen ("Messages in.txt", "w");
    if (F2==NULL) printf ("Messages in - NULL"), exit(2);

    int vector[24], message_in[12];

    for(int i=0;i<8288280;i++){
        for(int j=0;j<24;j++) fscanf(F1,"%d ", &(vector[j]));

        Decoding(vector,message_in);

        for(int j=0;j<12;j++) fprintf(F2,"%d ", (message_in[j]));
        fprintf(F2,"\n");
    }

    fclose(F1);
    fclose(F2);

    return 0;
}
