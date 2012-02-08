void populationCalculation(int maxSeqLen, int numEM, 
               Fitness *fitness, 
               int startPWMfound, 
               int minminSites, 
               double maxpFactor, 
               int numSeq, int numSeqEM, 
               char **seq, char **rseq, 
               int *seqLen, char *Iseq, 
               double *bfreq, 
               double **posWeight, int weightType, 
               double pvalueCutoff, int *emSeqLen, 
               double **pwm, int pwmLen, 
               double **epwm, double **opwm,
               char *pwmConsensus, int *scoreCutoff,
               char *sdyad, 
               int ii)
{


  int jj=0,llrDim;

  // two pwms before and after EM steps
  double **t1pwm=alloc_double_double(MAX_PWM_LENGTH,4);
  double **t2pwm=alloc_double_double(MAX_PWM_LENGTH,4);

  //These need to be allocated when running in parallel
  double **logpwm;              // log-transformed EM-derived EM
  double ** score, ** rscore; // subsequence score, plus and minus strands
  int **ipwm;      // integer pwm for computing llr score distribution
  Sites *siteEM;                // binding sites in EM sequences
  Pgfs *llrDist;                // llr distribution from pgf method

  llrDist=alloc_distr(MAX_DIMENSION);
  ipwm=alloc_int_int(MAX_PWM_LENGTH,4);
  logpwm=alloc_double_double(MAX_PWM_LENGTH,4);
  score=alloc_double_double(numSeq,maxSeqLen);
  rscore=alloc_double_double(numSeq,maxSeqLen);
  siteEM=alloc_site(MAX_SITES);

	
  int maxp=0, nsitesEM=0;
  double pwmDiff=0;

  if (!startPWMfound)
  {
    // This function modifles sdyad[ii]
    pwm_profile(pwm, pwmLen, sdyad);
  }
  // Make a copy and then subject the copy to EM
  copy_pwm(pwm, t1pwm, pwmLen);

  // standarize pwm
  standardize_pwm(t1pwm, pwmLen);
  // EM on randomly selected sequences
  maxp=(int)(maxpFactor*numSeqEM);

  //Check if the user wants to interrupt
  for(jj=0; jj<numEM; jj++)
  {
    // Modifies logpwm
    log_pwm(t1pwm,logpwm, pwmLen);
    // Compute ll score of each w-mer | motif model
    // Modifies score and rscore
    ll_score_motif_model(numSeq,seq,rseq,seqLen,logpwm,pwmLen,score,rscore,Iseq,bfreq);

    // compute p(zij|y=1) probability of binding sites started at position j on seq i
    // modifies score and rscore
    normalize(score,rscore,seqLen,pwmLen,numSeq,Iseq,maxp,posWeight,weightType);

    // E-step
    // Modifies t2pwm
    construct_pwm(t2pwm,score,rscore,seq,rseq,seqLen,numSeq,pwmLen,Iseq);

    // M-step
    // Modifies t2pwm
    standardize_pwm(t2pwm,pwmLen);
    // Only compute pwmDiff
    pwmDiff=check_convergence(t1pwm,t2pwm,pwmLen);
    // copy t2pwm to t1pwm
    copy_pwm(t2pwm,t1pwm,pwmLen);
    if (pwmDiff<=PWM_CONVERGENCE)  break;
  }

  // Copy t1pwm to epwm[ii]
  copy_pwm(t1pwm,epwm,pwmLen);
  // Modifies ipwm
  log_ratio_to_int(epwm,ipwm,pwmLen,bfreq);
  // Compute score distribution of the (int)PWM using Staden's method
  // Modifies llrDist
  llrDim=pwm_score_dist(ipwm,pwmLen,llrDist,bfreq);
  // Compute the score cutoff
  // Does not modify anything
  *scoreCutoff=determine_cutoff(llrDist,llrDim,pvalueCutoff);

  // test each w-mer to see if a motif site
  // test statistic: ll 
  // distribution: Staden method 
  // cutoff: user-specified
  // modifies siteEM
  nsitesEM=scan_em_seq_ptable(llrDist,llrDim,siteEM,numSeq,seq,rseq,seqLen,ipwm,pwmLen,*scoreCutoff,bfreq,Iseq);
   // loose threshould at this step, as em only on a subset of sequences
  if (nsitesEM>=max(2,minminSites))
  { 
    // construct pwm from the identified sites
    // Modifies opwm
    align_sites_count(siteEM,seq,rseq,nsitesEM,pwmLen,opwm);
    // Modifies opwm
    standardize_pwm(opwm,pwmLen);
    // Modifies pwmConsensus
    consensus_pwm(opwm,pwmLen,pwmConsensus);
    // compute E-value of the relative entroy score of each motif, use it as fitness
    // Only compute fitness.value
    (*fitness).value=E_value(opwm,nsitesEM,bfreq,pwmLen,numSeqEM,emSeqLen);
  }
  else
  {
    // if too few sites in a motif
    // Modifies opwm
    align_sites_count(siteEM,seq,rseq,nsitesEM,pwmLen,opwm);
    // Modifies opwm
    standardize_pwm(opwm,pwmLen);
    // Modifies pwmConsensus
    consensus_pwm(opwm,pwmLen,pwmConsensus);
    // Only compute fitness.value
    (*fitness).value=DUMMY_FITNESS;
  }
  (*fitness).index=ii;


  if (llrDist)
  { 
    free(llrDist);
    llrDist=NULL;
  }  
  if (siteEM)
  { 
    free(siteEM);
    siteEM=NULL;
  }  
  if (t1pwm[0])
  {
    free(t1pwm[0]);
    t1pwm[0]=NULL;
  }
  if (t1pwm)
  {
    free(t1pwm);
    t1pwm=NULL;
  }
  if (t2pwm[0])
  { free(t2pwm[0]);
    t2pwm[0]=NULL;
  }
  if (t2pwm)  
  { free(t2pwm);
    t2pwm=NULL;
  }
  if (logpwm[0])
  { free(logpwm[0]);
    logpwm[0]=NULL;
  }
  if (logpwm)
  { 
    free(logpwm);
    logpwm=NULL;
  }
  if(ipwm[0])
  {
    free(ipwm[0]);
    ipwm[0]=NULL;
  }
  if (ipwm)
  {
    free(ipwm);
    ipwm=NULL;
  }
  if (score[0])
  {
    free(score[0]);
    score[0]=NULL;
  }
  if (score)
  { 
    free(score);
    score=NULL;
  }
  if (rscore[0])
  { 
    free(rscore[0]);
    rscore[0]=NULL;
  }
  if (rscore)
  { 
    free(rscore);
    rscore=NULL;
  }
}
