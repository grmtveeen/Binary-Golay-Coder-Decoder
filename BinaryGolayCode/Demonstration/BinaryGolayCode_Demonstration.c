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

//Generate matrix G[B||I]
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

void Encoding (int u[], int result[]){   //"coder" gets message (mass) of 12 bits and converts they in vector (mass) of 24 bits
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
     int message_out[12]={1,0,0,1,1,1,0,1,0,1,1,0}, vector[24], message_in[12];

    //Message encoding and going to the channel
    Encoding(message_out,vector);

    printf("-----------------------------------------------------------------\n");
    printf("The message out: ");
    for(int i=0;i<12;i++) printf("%d ", message_out[i]);
    printf("\n-----------------------------------------------------------------\n");
    printf("Out vector:      ");
    for(int i=0;i<24;i++) printf("%d ", vector[i]);

    //Vector have 3 errors in channel
    if(vector[5]) vector[5]=0;
    else vector[5]=1;
    if(vector[13]) vector[13]=0;
    else vector[13]=1;
    if(vector[22]) vector[22]=0;
    else vector[22]=1;

    //Transmitted vector
    printf("\n-----------------------------------------------------------------\n");
    printf("In  vector:      ");
    for(int i=0;i<24;i++) printf("%d ", vector[i]);

    //Vector decodings
    Decoding(vector,message_in);

    //Get message
    printf("\n-----------------------------------------------------------------\n");
    printf("The message in:  ");
    for(int i=0;i<12;i++) printf("%d ", message_in[i]);
    printf("\n-----------------------------------------------------------------\n");

    return 0;
}
