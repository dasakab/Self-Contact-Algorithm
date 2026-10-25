/* Experimental linear-algebra-only acceleration. Not a contact algorithm change. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <superlu/slu_ddefs.h>

static SuperMatrix L,U;
static int *perm_r=NULL, *perm_c=NULL;
static int cached_n=0, cached=0;

static void forget(void){
  if(cached){ Destroy_SuperNode_Matrix(&L); Destroy_CompCol_Matrix(&U); }
  free(perm_r);free(perm_c);perm_r=NULL;perm_c=NULL;cached_n=0;cached=0;
}

void colsol_(double *a,int *nwk,int *maxa,int *nnm,double *v,int *nn,int *sflag,int *iflag){
  const int n=*nn;
  *iflag=0;
  if(*sflag==1){
    forget();
    size_t *count=(size_t*)calloc(n,sizeof(size_t));
    size_t nnz=0;
    if(!count){fprintf(stderr,"SPARSE: count alloc error\n");exit(90);}
    for(int i=0;i<n;i++){
      int pos=maxa[i]-1;
      int height=maxa[i+1]-maxa[i];
      if(pos<0 || pos+height>*nwk){fprintf(stderr,"SPARSE: invalid skyline index i=%d pos=%d height=%d nwk=%d\n",i,pos,height,*nwk);exit(91);}
      for(int k=0;k<height;k++){
        double x=a[pos+k];
        if(x!=0){ int j=i-k; if(j<0){fprintf(stderr,"SPARSE: skyline height invalid\n");exit(92);} count[j]++; if(j!=i) count[i]++; nnz+=(j!=i)?2:1; }
      }
    }
    int_t *colptr=(int_t*)malloc((n+1)*sizeof(int_t));
    int_t *rows=(int_t*)malloc(nnz*sizeof(int_t));
    double *values=(double*)malloc(nnz*sizeof(double));
    int_t *next=(int_t*)malloc(n*sizeof(int_t));
    if(!colptr||!rows||!values||!next){fprintf(stderr,"SPARSE: matrix alloc error\n");exit(93);}
    colptr[0]=0;
    for(int i=0;i<n;i++)colptr[i+1]=colptr[i]+count[i];
    memcpy(next,colptr,n*sizeof(int_t));
    for(int i=0;i<n;i++){
      int pos=maxa[i]-1;
      int height=maxa[i+1]-maxa[i];
      for(int k=0;k<height;k++){
        double x=a[pos+k];
        if(x==0)continue;
        int j=i-k; int_t at=next[j]++;
        rows[at]=i; values[at]=x;
        if(j!=i){at=next[i]++;rows[at]=j;values[at]=x;}
      }
    }
    free(next);free(count);
    fprintf(stderr,"SPARSE_COLSOL: n=%d, skyline=%d, sparse entries=%zu, factorizing\n",n,*nwk,nnz);fflush(stderr);
    SuperMatrix A,B;
    double *zero=(double*)calloc(n,sizeof(double));
    if(!zero){fprintf(stderr,"SPARSE: zero alloc error\n");exit(94);}
    dCreate_CompCol_Matrix(&A,n,n,(int_t)nnz,values,rows,colptr,SLU_NC,SLU_D,SLU_GE);
    dCreate_Dense_Matrix(&B,n,1,zero,n,SLU_DN,SLU_D,SLU_GE);
    perm_r=(int*)malloc(n*sizeof(int));perm_c=(int*)malloc(n*sizeof(int));
    superlu_options_t options;
    set_default_options(&options);
    options.ColPerm=COLAMD;
    SuperLUStat_t stat;StatInit(&stat);
    int_t info=0;
    dgssv(&options,&A,perm_c,perm_r,&L,&U,&B,&stat,&info);
    StatFree(&stat);
    Destroy_SuperMatrix_Store(&A);
    Destroy_SuperMatrix_Store(&B);
    free(rows); free(colptr); free(values);free(zero);
    if(info!=0){fprintf(stderr,"SPARSE_COLSOL: factorisation FAIL info=%ld\n",(long)info);*iflag=1;forget();return;}
    cached_n=n;cached=1;
    fprintf(stderr,"SPARSE_COLSOL: factorisation complete\n");fflush(stderr);
  }else if(*sflag==2){
    if(!cached||n!=cached_n){fprintf(stderr,"SPARSE_COLSOL: no saved factorisation!\n");exit(95);}
    SuperMatrix B;dCreate_Dense_Matrix(&B,n,1,v,n,SLU_DN,SLU_D,SLU_GE);
    SuperLUStat_t stat;StatInit(&stat);int info=0;
    dgstrs(NOTRANS,&L,&U,perm_c,perm_r,&B,&stat,&info);
    StatFree(&stat);Destroy_SuperMatrix_Store(&B);
    if(info!=0){fprintf(stderr,"SPARSE_COLSOL: solve FAIL info=%d\n",info);*iflag=1;}
  }else{ fprintf(stderr,"SPARSE_COLSOL: unknown flag %d\n",*sflag); exit(96); }
}
