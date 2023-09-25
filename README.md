# dynamoDMT v0.2b

Filament analysis using Dynamo

- aa_setup.m (CPU)

- aa_correctFlaDirect.m (CPU) * Important for axoneme not other

- aa_imodModel2Filament.m (CPU)

- aa_cropAndAverage.m (CPU)

- aa_intraAln.m (GPU)

- aa_alignIntraAvg.m (CPU)

- aa_alignAllParticles.m (GPU)

- aa_generateAxonemeAverage.m (CPU) * Optional for axoneme

- aa_filamentRepick.m (CPU)

- aa_intraAvgRepick.m (GPU)

- aa_alignRepickAvg.m (CPU)

- aa_alignRepickParticles (GPU)



For doublet MTD, make sense to first align with 8-nm repeat, then repick with 16-nm repeat
