On the slides: 

- Add bibliography w/ paper screenshots
- Mention GOF tests w/ U-CARE and R2ucare, chapter on Gentle introduction
- Add some of our slides on multievent models
- Add a slide Summary and/or Take-home messages
- Add a slide w/ bibliography and/or To go further
- Write a slide w/ tricks to speed up the computations and improve convergence
    - Marginalization (cite relevant paper I review for Ecol Monog); by hand, or using what is done in Nimble/Stan. 
    - Nimble w/ NimbleEcology, Stan w/ tag Ecology
    - Provide examples of IPM w/ Stan, with Nimble, with TMB/ADMB 
    - m-arrays, but limitations when it comes to individual effect (random of covariates)
    - TMB, Rcpp; to get you started, I have some code
    - Other stuff?
- Introduce HMMs and our Ecol Let papers; actually, we use HMMs not general SSMs for IPMs; see also recent paper by Takis in Biometrics formulating IPMs as HMMs!
- Present CJS as a particular case of the multistate, as in https://github.com/oliviergimenez/multievent_jags_R
- Illustrate what dcat does
- Add a few memes

On the practicals:

- Add an example of multievent model in the multistate R script? Use the roe deer example from my TPB paper; check out the examples we cover in 
- Put m-array stuff in another file, gathering single and multistate stuff (or two files, yet to be decided)
- Provide link to Nimble translation of book? Add Nimble implementation of the exercises?
- Translate the R scripts in Rmd?
- The multistate example can be simplified, we don't need the loops on individuals/time to define the state matrices. In the exercise, we need the time loop. Mention that to add individual stuff, an extra loop on indiviuals is needed. Tension is with the lecture that presents the model with time and individual effect as the most general and generic model
- Add multievent model https://github.com/oliviergimenez/multievent_jags_R ; make it an exercise; 
