//=======================================================================
// This is part of the 2DECOMP&FFT library
// 
// 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
// decomposition. It also implements a highly scalable distributed
// three-dimensional Fast Fourier Transform (FFT).
//
// Copyright (C) 2009-2021 Ning Li, the Numerical Algorithms Group (NAG)
//
//=======================================================================

// This is the shared-memory code using System V IPC API

/*
 This shared-memory code is kindly provided by David Tanqueray of Cray Inc.
 who also helped the author adapt it to use in 2DECOMP&FFT. His assistance
 is greatly appreciated.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <mpi.h>

#ifndef DBG
#define DBG
#endif
static int shm_debug=0;

float log2f(float); 

void set_shm_debug_();
void get_smp_map2_(MPI_Fint *comm, MPI_Fint *nnodes, MPI_Fint *my_node,
                  MPI_Fint *ncores, MPI_Fint *my_core, MPI_Fint *maxcor);
void alloc_shm_(MPI_Aint *ptr, MPI_Fint *nelem, MPI_Fint *type,
                MPI_Fint *comm, MPI_Fint *ret);
void dealloc_shm_(MPI_Aint *ptr, MPI_Fint *comm);

void set_shm_debug_()
{
    shm_debug = 1;
}

void get_smp_map2_(MPI_Fint *comm, MPI_Fint *nnodes, MPI_Fint *my_node,
                  MPI_Fint *ncores, MPI_Fint *my_core, MPI_Fint *maxcor)
{
    MPI_Comm world;
    int err, pe, mype, npes, nnds, ncrs, maxc;
    int nlen, mynid, *nidlst, *nodlst;
    int i, n;
    char nodnam[MPI_MAX_PROCESSOR_NAME];
    char string[20];
    FILE *fp;

    MPI_Comm_rank(MPI_COMM_WORLD,&pe);

    world = MPI_Comm_f2c(*comm);
    MPI_Comm_rank(world, &mype);
    MPI_Comm_size(world, &npes);
    MPI_Get_processor_name(nodnam, &nlen);
#ifdef USE_NAME
    mynid = atoi(nodnam+3);
#else
    sprintf(string," pe %d /proc/cray_xt/nid",mype);
    if ((fp = fopen("/proc/cray_xt/nid", "r")) ==  NULL) {
          perror(string);
          exit(1);
    }
    fscanf(fp,"%i", &mynid);
    fclose(fp);
#endif
#ifdef DBG
    if (shm_debug) {
      fprintf(stderr," pe %d mype %d of %d, nodnam = %s (len = %d), node = %d\n",
             pe, mype, npes, nodnam, nlen, mynid);
      MPI_Barrier(world);
    }
#endif

    /* get list of nodeid for each pe */
    nidlst = malloc(npes*sizeof(int));
    MPI_Allgather(&mynid, 1, MPI_INT, nidlst, 1, MPI_INT, world);

    nodlst = malloc(npes*sizeof(int));
    nnds = ncrs = 0;
    for (i=0; i<npes; i++) {
      /* get my core id and no. of cores on my node */
      if (i == mype) *my_core = (MPI_Fint)(ncrs);
      if (nidlst[i] == mynid) ncrs++;

      /* get my node id and number of unique nodes */
      for (n=0; n<nnds; n++) {
        if (nodlst[n] == nidlst[i]) break;            /* node already in list */
      }
      if (nidlst[i] == mynid) *my_node = (MPI_Fint)(n); /* save my node index */
      if (n >= nnds) nodlst[nnds++] = nidlst[i];      /* add new node to list */
    }
    /* get max core counts over all nodes */
    MPI_Allreduce(&ncrs, &maxc, 1, MPI_INT, MPI_MAX, world);
    *nnodes = (MPI_Fint)(nnds);
    *ncores = (MPI_Fint)(ncrs);
    *maxcor = (MPI_Fint)(maxc);

#ifdef DBG
    if (shm_debug) {
      fprintf(stderr," pe %d nnodes=%d ncores=%d maxcor=%d\n",pe,nnds,ncrs,maxc);
      for (n=0; n<nnds; n++) {
        if (nodlst[n]==mynid && *my_core==0) {
          fprintf(stderr," pe %d node %d (%s) of %d nodes, with %d cores\n",
                 pe,n,nodnam,nnds,ncrs);
        }
        MPI_Barrier(world);
      }
      fprintf(stderr,"max core count = %d\n",maxc);
    }
#endif
    free(nidlst);
    free(nodlst);
}


/* Implementation using System V shared memory [shmget & shmat] */
/* Because this uses 32-bit sizes, we allow for multiple segments */

#define BLKSIZE (1L<<30)

struct blktyp {
    struct blktyp *next;
    size_t size;
    void *addr;
    int shmid;
    key_t key;
};
typedef struct blktyp Blktyp;

struct segtyp {
    struct segtyp *next;
    void *base;
    Blktyp *blks;
};
typedef struct segtyp Segtyp;

Segtyp *seglst=NULL;

void alloc_shm_(MPI_Aint *ptr, MPI_Fint *nelem, MPI_Fint *type,
                MPI_Fint *comm, MPI_Fint *ret)
{
    MPI_Comm world;
    MPI_Fint err;
    int pe, np, mype, typsiz, msk;
    size_t size, blksize;
    char *shm, *shmaddr;
    char string[20];
    int shmid;
    key_t key;
    Segtyp *seg;
    Blktyp *blk;
#ifdef DBG
    int npes,*plst,i;
    char plist[1024],*s;
#endif

    MPI_Comm_rank(MPI_COMM_WORLD, &pe);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    msk = ( 1 << ((int)ceilf(log2f((float)np))) ) - 1;

    world = MPI_Comm_f2c(*comm);
    MPI_Comm_rank(world, &mype);
#ifdef DBG
    if (shm_debug) {
      MPI_Comm_size(world, &npes);
      plst = (int*)malloc(npes*sizeof(int));
      MPI_Allgather(&pe,1,MPI_INT,plst,1,MPI_INT,world);
      for (i=0,s=plist; i<npes; i++) s=s+sprintf(s,"%d,",plst[i]);
      s[0]='\0'; if (s>plist) s[-1]='\0';
      fprintf(stderr," pe %d al_shm: comm sz/rk=%d/%d pes=%s\n",pe,npes,mype,plist);
    }
#endif
    MPI_Type_size(MPI_Type_f2c(*type), &typsiz);
    size = (size_t)*nelem * (size_t)typsiz;

    err = 0;
    /* setup structure to keep track of this segment and its blocks */
    seg = (Segtyp *)malloc(sizeof(Segtyp));
    seg->next = seglst;
    seg->base = NULL;
    seg->blks = NULL;
    seglst = seg;

    shmaddr = NULL;
    while (size) {
      blksize = size<BLKSIZE ? size : BLKSIZE;

      if (!mype) key = rand()&~msk | pe&msk;
      MPI_Bcast(&key, 1, MPI_INT, 0, world);

      sprintf(string," pe %d shmget",pe);
      if ((shmid = shmget(key, blksize, IPC_CREAT | 0666)) < 0) {
        perror(string);
        if (*ret) err++; else exit(1);
      }
#ifdef DBG
      if (shm_debug) fprintf(stderr," pe %d shmget: key = %d, shmid = %d\n",
                            pe,key,shmid);
#endif

      sprintf(string," pe %d shmat",pe);
      if ((shm = shmat(shmid, shmaddr, 0)) == (char *) -1) {
        perror(string);
        err++;
      }

      /* wait for all cores to attach block before ... */
      MPI_Barrier(world);

      /* ... marking it for deletion */
      sprintf(string," pe %d shmctl",pe);
      if (shmctl(shmid, IPC_RMID, NULL) < 0) {
        perror(string);
        if (*ret) err++; else exit(1);
      }
#ifdef DBG
      if (shm_debug) fprintf(stderr," pe %d shmat: shm = %lx, size = %d, limit = %lx\n",pe,shm,blksize,shm+blksize);
#endif
      if (seg->base == NULL) seg->base = shm;

      /* setup structure to record block details */
      blk = malloc(sizeof(Blktyp));
      blk->size = blksize;
      blk->addr = shm;
      blk->shmid = shmid;
      blk->key   = key;
      /* and add it to segment's list */
      blk->next = seg->blks;
      seg->blks = blk;

      shmaddr = shm + blksize;
      size -= blksize;
    }

    *ptr = (MPI_Aint)seg->base;
    if (*ret) *ret = err;
    else if (err) exit(1);
}

void dealloc_shm_(MPI_Aint *ptr, MPI_Fint *comm)
{
    MPI_Comm world;
    int pe, mype, err;
    char string[20];
    Segtyp *seg;
    Segtyp **prev;
    Blktyp *blk=NULL, *nblk;
    void *shm;

    MPI_Comm_rank(MPI_COMM_WORLD, &pe);
#ifdef DBG
    if (shm_debug) fprintf(stderr," pe %d dealloc_shm: ptr=%lx\n",pe,*ptr);
#endif

    world = MPI_Comm_f2c(*comm);
    MPI_Comm_rank(world, &mype);
    err = 0;
    /* Find segment with specified start address and remove from list */
    seg = seglst;
    prev = &seglst;
    while (seg) {
      if (seg->base == (void*)*ptr) {
        blk=seg->blks;
        *prev = seg->next;
        free(seg);
        break;
      }
      prev = &seg->next;
      seg = seg->next;
    }
    if (blk == NULL) {
      fprintf(stderr," pe %d dealloc_shm: segment at address %lx not found\n",
              pe, *ptr);
    }
    /* detach all blocks in this segment */
    while (blk) {
      shm = blk->addr;
      sprintf(string," pe %d shmdt",pe);
      if (shmdt((char*)shm) < 0) {
        perror(string);
        err++;
      }
      nblk = blk->next;
      free(blk);
      blk = nblk;
    }
    if (err) exit(1);
}

