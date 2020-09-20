#include<bits/stdc++.h> 

void getCI(vector<double> &list);
void hits();

double pav, CImin, CImax, sustop, susbot, inftop, infbot;
double sustopav = 0, susbotav = 0, inftopav = 0, infbotav = 0;
double sustopav2 = 0, susbotav2 = 0, inftopav2 = 0, infbotav2 = 0;

void doreps(string file)
{
  const short varmax = 16;
  
  const long  nrep = 40;
  long rep, rep2, i, v;
  long nsire, M;
  double avq, avq2, avqbre, avqbre2, avqqbre, te;
  double av = 0, av2 = 0, avb = 0, avb2 = 0, avc = 0, avc2 = 0, ac, acb, acc;
  vector< vector <double> > varst;
  
  long timest;
  double varav[varmax][nrep], varCImin[varmax][nrep], varCImax[varmax][nrep];
  
  if(1==1){  // This gets results for the first SIRE paper
		const int vmax = 15;
		double val;
		double vmean[vmax], vmean2[vmax],vmeanav[vmax],vmeanav2[vmax], vsdav[vmax], vsdav2[vmax];
		string name, vname[vmax];
		
		for(v = 0; v < vmax; v++){ vmeanav[v] = 0; vmeanav2[v] = 0; vsdav[v] = 0;  vsdav2[v] = 0;}
			
		for(rep = 0; rep < nrep; rep++){
			cout << "rep: "<< rep << "\n";
			stringstream ss;
			ss << file << "_" << rep << ".xml";
			cout << ss.str() << "\n";

			readfile(ss.str());                                             // Reads in the input file

			init(); 
	
			for(v = 0; v < vmax; v++){ vmean[v] = 0; vmean2[v] = 0;}
			for(s = 0; s < nsamp; s++){
				//cout << s << " samp\n";
				update();                                // Performs MCMC sampling
				
				if(s >= burnin){
					for(v = 0; v < vmax; v++){
						switch(v){
						case 0: val = beta; name="beta"; break;
						case 1: val = gama; name="gamma"; break;
						case 2: val = kshape; name="k";  break;
						case 3: val = a_g; name="a_g";  break;
						case 4: val = a_f; name="a_f";  break;
						case 5: val = a_r; name="a_r";  break;
						case 6: val = delta_g; name="delat_g";  break;
						case 7: val = delta_f; name="delat_f";  break;
						case 8: val = delta_r; name="delat_r";  break;
						case 9: val = vare_gg; name="vare_gg";  break;
						case 10: val = vare_ff; name="vare_gg";  break;
						case 11: val = vare_rr; name="vare_gg";  break;
						case 12: val = vare_gf; name="vare_gg";  break;
						case 13: val = vare_gr; name="vare_gg";  break;
						case 14: val = vare_fr; name="vare_gg";  break;
						}
						vname[v]=name;
						vmean[v] += val; vmean2[v] += val*val; 
					}
				}
			}
			
			for(v = 0; v < vmax; v++){
				vmean[v] /= (nsamp-burnin);
				vmean2[v] /= (nsamp-burnin);
				vmeanav[v] += vmean[v];
				vmeanav2[v] += vmean[v]*vmean[v];
				
				double sd =  sqrt(vmean2[v]-vmean[v]*vmean[v]);
				vsdav[v] += sd;
				vsdav2[v] += sd*sd;
			}
		
			for(v = 0; v < vmax; v++){
				double mean = vsdav[v]/(rep+1);
				double dif = sqrt(vsdav2[v]/(rep+1) -mean*mean);
				cout << vname[v] << " "<< vmeanav[v]/(rep+1) << " " << mean << " " << mean-dif << " " << mean+dif << "\n";
			}
			
			for(v = 3; v < 6; v++){
				double mean = vmeanav[v]/(rep+1);
				double dif = sqrt( vmeanav2[v]/(rep+1) -mean*mean);
				cout << mean << " " << mean-dif << " " << mean+dif << " ";
			}
			cout << "  means\n";
			
			for(v = 3; v < 6; v++){
				double mean = vsdav[v]/(rep+1);
				double dif = sqrt(vsdav2[v]/(rep+1) -mean*mean);
				cout << mean << " " << mean-dif << " " << mean+dif << " ";
			}
			cout << "  SD\n";
		}
		
		return;	
  }
  
  timest = clock();
  for(rep = 0; rep < nrep; rep++){
    cout << "rep: "<< rep << "\n";
    stringstream ss;
	ss << file << "_" << rep << ".xml";
	//ss <<  "Simulate/SIM/sim1337_19.xml";
    cout << ss.str() << "\n";

    readfile(ss.str());                                             // Reads in the input file

    init();                                                     // Initialises MCMC quantites

	/*
	vector <long> sirelist;
	if(sirelistst.size() > 0){
		long j;
	
		for(i = 0; i < sirelistst.size(); i++){
			j = long(atof( sirelistst[i].c_str()));
			 sirelist.push_back(j);
		}
	}
	else{
	    nsire = 0; while(nsire < N && indtrial[nsire] == -1) nsire++;
		for(i = 0; i < nsire; i++) sirelist.push_back(i);
		//for(i = nsire; i < N; i++) sirelist.push_back(i);
	}
	*/
	
    //susp.clear(); infp.clear(); sush2.clear(); infh2.clear(); betast.clear();
    varst.clear(); varst.resize(varmax);
  
	cout << beta << " beta \n";
	cout << gama << " gam \n";
	cout << gama << " kshape \n";
	return;
   
	
    timetot -= clock();
    for(s = 0; s < nsamp; s++){
		cout << s << " samp\n";
       update();                                // Performs MCMC sampling
      
      if(s >= burnin){
	varst[0].push_back(vara_gg+vare_gg);
	varst[1].push_back(vara_ff+vare_ff);
	varst[2].push_back(vara_rr+vare_rr);
	varst[3].push_back(vara_gg/(vara_gg+vare_gg));
	varst[4].push_back(vara_ff/(vara_ff+vare_ff));
	varst[5].push_back(vara_rr/(vara_rr+vare_rr));
	varst[6].push_back(vara_gf/sqrt(vara_gg*vara_ff));
	varst[7].push_back(vare_gf/sqrt(vare_gg*vare_ff));
	varst[8].push_back(vara_gr/sqrt(vara_gg*vara_rr));
	varst[9].push_back(vare_gr/sqrt(vare_gg*vare_rr));
	varst[10].push_back(vara_fr/sqrt(vara_rr*vara_ff));
	varst[11].push_back(vare_fr/sqrt(vare_rr*vare_ff));
	varst[12].push_back(beta);
	varst[13].push_back(gama);
	varst[14].push_back(kshape);
	 varst[15].push_back(siggeff_g);
	
      }
      
      if(s > burnin){
	nqav++; for(i = 0; i < N; i++){ q_g_av[i] += q_g[i]; q_f_av[i] += q_f[i]; if(mod == SIR) q_r_av[i] += q_r[i];}
      }
      
      if(s > burnin && (s%1000 == 0 || s == nsamp-1)){
	//cout << rep << " " << s << "\n";
	/*
	for(loop = 1; loop >= 0; loop--){
	  switch(loop){
	    case 0: imin = 0; imax = nsire; break;
	    case 1: imin = nsire; imax = N; break;
	  }
	   switch(loop){
	    case 0: cout << "results for sires: ";  break;
	    case 1: cout << "results for prog: ";  break;
	  }
	  */
	  /*
	  M =sirelist.size();
	long k;
	  avq = 0; avq2 = 0; avqbre = 0; avqbre2 = 0; avqqbre = 0;
	  for(k = 0; k < M; k++){
		  i= sirelist[k];
	    avq += q_g_av[i]/nqav; avq2 += (q_g_av[i]/nqav)*(q_g_av[i]/nqav); 
	    avqbre += q_g_bv[i]; avqbre2 += q_g_bv[i]*q_g_bv[i]; 
	    avqqbre += (q_g_av[i]/nqav)*q_g_bv[i];
	  }
	  ac = (avqqbre/M - (avq/M)*(avqbre/M))/(sqrt(avq2/M - (avq/M)*(avq/M))*sqrt(avqbre2/M - (avqbre/M)*(avqbre/M)));
	  cout << ac << " ";

	  avq = 0; avq2 = 0; avqbre = 0; avqbre2 = 0; avqqbre = 0;
	  for(k = 0; k < M; k++){
		  i= sirelist[k];
	    avq += q_f_av[i]/nqav; avq2 += (q_f_av[i]/nqav)*(q_f_av[i]/nqav); 
	    avqbre += q_f_bv[i]; avqbre2 += q_f_bv[i]*q_f_bv[i]; 
	    avqqbre += (q_f_av[i]/nqav)*q_f_bv[i];
	  }
	  acb = (avqqbre/M - (avq/M)*(avqbre/M))/(sqrt(avq2/M - (avq/M)*(avq/M))*sqrt(avqbre2/M - (avqbre/M)*(avqbre/M)));
	 
	  cout << acb << " ";
	  
	  if(q_r_bv.size() > 0){
	    avq = 0; avq2 = 0; avqbre = 0; avqbre2 = 0; avqqbre = 0;
	     for(k = 0; k < M; k++){
		  i= sirelist[k];
	      avq += q_r_av[i]/nqav; avq2 += (q_r_av[i]/nqav)*(q_r_av[i]/nqav); 
	      avqbre += q_r_bv[i]; avqbre2 += q_r_bv[i]*q_r_bv[i]; 
	      avqqbre += (q_r_av[i]/nqav)*q_r_bv[i];
	    }
	    acc = (avqqbre/M - (avq/M)*(avqbre/M))/(sqrt(avq2/M - (avq/M)*(avq/M))*sqrt(avqbre2/M - (avqbre/M)*(avqbre/M)));
	  }
	  
	  cout << acc << "  a\n";
	  */
      }

    }
	
		
		/*
    ofstream corplot("Simulate/GP/cor3"); 
    for(i = 0; i < nsire; i++){
      corplot << q_g_bv[i] << " " << q_g_av[i]/nqav << " "<< q_f_bv[i] << " " << q_f_av[i]/nqav << " "<< q_r_bv[i] << " " << q_r_av[i]/nqav << "\n";
    }
   
   
    ofstream corprogplot("Simulate/GP/corprog3"); 
    for(i = nsire; i < N; i++){
      corprogplot << q_g_bv[i] << " " << q_g_av[i]/nqav << " "<< q_f_bv[i] << " " << q_f_av[i]/nqav << " "<< q_r_bv[i] << " " << q_r_av[i]/nqav << "\n";
    }
    */
	
   // hits();
      
    cout << rep << "rep\n";
    cout << ac << " "<< acb << " " << acc << " acc\n";
    
    for(v = 0; v < varmax; v++){
		
		varst[0].push_back(vara_gg+vare_gg);
	varst[1].push_back(vara_ff+vare_ff);
	varst[2].push_back(vara_rr+vare_rr);
	varst[3].push_back(vara_gg/(vara_gg+vare_gg));
	varst[4].push_back(vara_ff/(vara_ff+vare_ff));
	varst[5].push_back(vara_rr/(vara_rr+vare_rr));
	varst[6].push_back(vara_gf/sqrt(vara_gg*vara_ff));
	varst[7].push_back(vare_gf/sqrt(vare_gg*vare_ff));
	varst[8].push_back(vara_gr/sqrt(vara_gg*vara_rr));
	varst[9].push_back(vare_gr/sqrt(vare_gg*vare_rr));
	varst[10].push_back(vara_fr/sqrt(vara_rr*vara_ff));
	varst[11].push_back(vare_fr/sqrt(vare_rr*vare_ff));
	varst[12].push_back(beta);
	varst[13].push_back(gama);
	varst[14].push_back(kshape);
	 varst[15].push_back(siggeff_g);
	 
      switch(v){
	case 0: cout << "sussize"; break;
	case 1: cout << "infsize"; break;
	case 2: cout << "sush2"; break;
	case 3: cout << "infh2"; break;
	case 4: cout << "cora"; break;
	case 5: cout << "core"; break;
	case 6: cout << "beta"; break;
	case 7: cout << "susfi"; break;
	case 8: cout << "inffi"; break;
	case 9: cout << "sigma"; break;
      }
      cout << ":";
      getCI(varst[v]); varav[v][rep] = pav; varCImin[v][rep] = CImin; varCImax[v][rep] = CImax;
      cout << "\n";
    }
  
    av += ac; av2 += ac*ac;
    avb += acb; avb2 += acb*acb;
    avc += acc; avc2 += acc*acc; 
  }
  cout << "average correlations  sus/inf:  ";
  cout << av/nrep << " " <<  sqrt(av2/nrep-(av/nrep)*(av/nrep)) << " ";
  cout << avb/nrep << " " <<  sqrt(avb2/nrep-(avb/nrep)*(avb/nrep)) << " ";
  cout << avc/nrep << " " <<  sqrt(avc2/nrep-(avc/nrep)*(avc/nrep)) << " sus inf rec ans\n";
  
  for(v = 0; v < varmax; v++){
    for(rep = 0; rep < nrep; rep++){
      for(rep2 = 0; rep2 < nrep-1; rep2++){
	if(varav[v][rep2] > varav[v][rep2+1]){
	  te = varav[v][rep2]; varav[v][rep2] = varav[v][rep2+1]; varav[v][rep2+1] = te;
	  te = varCImin[v][rep2]; varCImin[v][rep2] = varCImin[v][rep2+1]; varCImin[v][rep2+1] = te;
	  te = varCImax[v][rep2]; varCImax[v][rep2] = varCImax[v][rep2+1]; varCImax[v][rep2+1] = te;
	}
      }
    }
  }
  
  /*
  for(rep = 0; rep < nrep; rep++){
    cout << rep << " ";
    for(v = 0; v < varmax; v++){
      cout << varav[v][rep] << " "<< varCImin[v][rep] << " "<< varCImax[v][rep] << " ";
    }
    cout << "\n";
  }
  */
  
   cout << "HITS result sustop susbot inftop infbot:  " 
	<< sustopav/nrep << " " << sqrt(sustopav2/nrep - (sustopav/nrep)*(sustopav/nrep)) << " "
	<< susbotav/nrep << " " << sqrt(susbotav2/nrep - (susbotav/nrep)*(susbotav/nrep)) << " "
	<< inftopav/nrep << " " << sqrt(inftopav2/nrep - (inftopav/nrep)*(inftopav/nrep)) << " "
	<< infbotav/nrep << " " << sqrt(infbotav2/nrep - (infbotav/nrep)*(infbotav/nrep)) << "\n";

  cout <<  (clock()-timest)/(CLOCKS_PER_SEC*60.0) << " time\n";
}

void getCI(vector<double> &list)
{
  long n, j, v; 
  double  f;
  
  n = list.size();
  if(n == 0){ pav = -1; CImin = -1; CImax = -1; return;}
  
  sort(list.begin(),list.end());
  pav = 0; for(j = 0; j < n; j++) pav += list[j]; pav /= n;

  v = long((n-1)*0.05); f = (n-1)*0.05 - v; CImin = list[v]*(1-f) +  list[v+1]*f;	
  v = long((n-1)*0.95); f = (n-1)*0.95 - v; CImax = list[v]*(1-f) +  list[v+1]*f;
 
  cout << pav << " " << CImin << " " << CImax << "  CI";
}

void hits()
{
  long i, j, ii;
  double te, num, numcor;
  long avlist[N], bvlist[N];

  for(i = 0; i < N; i++){ avlist[i] = i; bvlist[i] = i;}
    
  for(j = 0; j < N; j++){
    for(i = 0; i < N-1; i++){    // works out if top hits are correct
      if(q_g_av[i] > q_g_av[i+1]){
	te = q_g_av[i]; q_g_av[i] = q_g_av[i+1]; q_g_av[i+1] = te;
	ii = avlist[i]; avlist[i] = avlist[i+1]; avlist[i+1] = ii;
      }
      
      if(q_g_bv[i] > q_g_bv[i+1]){
	te = q_g_bv[i]; q_g_bv[i] = q_g_bv[i+1]; q_g_bv[i+1] = te;
	ii = bvlist[i]; bvlist[i] = bvlist[i+1]; bvlist[i+1] = ii;
      }
    }
  }
  
  num = 0; numcor = 0;
  for(i = long(0.9*N); i < N; i++){
    for(j = long(0.9*N); j < N; j++) if(bvlist[j] == avlist[i]){ numcor++; break;}
    num++; 
  }
  sustop = numcor/num;
  
  num = 0; numcor = 0;
  for(i = 0; i < long(0.1*N); i++){
    for(j = 0; j < long(0.1*N); j++) if(bvlist[j] == avlist[i]){ numcor++; break;}
    num++; 
  }
  susbot = numcor/num;
    
  
  // Infectivity 
  for(i = 0; i < N; i++){ avlist[i] = i; bvlist[i] = i;}
    
  for(j = 0; j < N; j++){
    for(i = 0; i < N-1; i++){    // works out if top hits are correct
      if(q_f_av[i] > q_f_av[i+1]){
	te = q_f_av[i]; q_f_av[i] = q_f_av[i+1]; q_f_av[i+1] = te;
	ii = avlist[i]; avlist[i] = avlist[i+1]; avlist[i+1] = ii;
      }
      
      if(q_f_bv[i] > q_f_bv[i+1]){
	te = q_f_bv[i]; q_f_bv[i] = q_f_bv[i+1]; q_f_bv[i+1] = te;
	ii = bvlist[i]; bvlist[i] = bvlist[i+1]; bvlist[i+1] = ii;
      }
    }
  }
  
  num = 0; numcor = 0;
  for(i = 0; i < long(0.1*N); i++){
    for(j = 0; j < long(0.1*N); j++) if(bvlist[j] == avlist[i]){ numcor++; break;}
    num++; 
  }
  infbot = numcor/num;
    
  num = 0; numcor = 0;
  for(i = long(0.9*N); i < N; i++){
    for(j = long(0.9*N); j < N; j++) if(bvlist[j] == avlist[i]){ numcor++; break;}
    num++; 
  }
  inftop = numcor/num;
  
 
  cout << "HITS sustop susbot inftop infbot:  " << sustop << " " << susbot << " " << inftop << " " << infbot << "p\n";
  sustopav += sustop; susbotav += susbot; inftopav += inftop; infbotav += infbot;
  sustopav2 += sustop*sustop; susbotav2 += susbot*susbot; inftopav2 += inftop*inftop; infbotav2 += infbot*infbot;
}
