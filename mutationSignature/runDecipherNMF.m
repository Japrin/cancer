function runDecipherNMF(workDir,totalSig)

  matInFile=strcat(workDir,'/SIn.mat')
  matOutFile=strcat(workDir,'/SOut.mat')
  
  cancerType='breast'
  sampleNames=textread(strcat(workDir,'/sampleNames.txt'),'%s','headerlines',1)
  types=textread(strcat(workDir,'/types.txt'),'%s','headerlines',1)
  subtypes=textread(strcat(workDir,'/subtypes.txt'),'%s','headerlines',1)
  originalGenomes=dlmread(strcat(workDir,'/originalGenomes.txt'))
  save(matInFile)
  
   
  addpath('/Share/BP/zhenglt/01.bin/mutation/WTSI_Mutational_Signature_Framework/source/');
  addpath('/Share/BP/zhenglt/01.bin/mutation/WTSI_Mutational_Signature_Framework/plotting/');
  
  %% Open matlabpool
  if ( matlabpool('size') == 0 )
      matlabpool open; % opens the default matlabpool, if it is not already opened
  end
  
  %% Define parameters 
  totalSignatures = totalSig;
  %%totalSignatures = 5;
  iterationsPerCore = 10;
  inputFile = matInFile;
  outputFile = matOutFile;
  
  %% Decipher the signatures of mutational processes from catalogues of mutations
  [input allProcesses allExposures idx processes exposures processStab processStabAvg] = ...
      decipherMutationalProcesses(iterationsPerCore, totalSignatures, inputFile, outputFile);
  
  %% Plotting the signatures of mutational processes
  plotSignatures(processes, input, allProcesses, idx, processStabAvg);
  %% Plotting the signature exposures
  plotSignaturesExposureInSamples(exposures, input);
  
  dlmwrite(strcat(workDir,'/S.H'),exposures,'\t');
  dlmwrite(strcat(workDir,'/S.W'),processes,'\t');

end
