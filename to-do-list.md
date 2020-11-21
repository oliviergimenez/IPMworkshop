On the slides: 

- Add some of our slides on multievent models
- Add a slide Summary and/or Take-home messages
- Write a slide w/ tricks to speed up the computations and improve convergence
    - Marginalization (cite relevant paper I review for Ecol Monog); by hand <https://esajournals.onlinelibrary.wiley.com/doi/10.1002/eap.2112>, or relying what is done in Nimble (see <https://cran.r-project.org/web/packages/nimbleEcology/index.html>) or Stan (see <https://mc-stan.org/docs/2_22/stan-users-guide/latent-discrete-chapter.html>). 
    - m-arrays, but limitations when it comes to individual effect (random of covariates)
    - Stan, Nimble are often faster than Jags, in that they need less MCMC samples to reach convergence; TMB, Rcpp; to get you started, I have some code
    - Other stuff?
- Introduce HMMs and our Ecol Let papers; actually, we use HMMs not general SSMs for IPMs; see also recent paper by Takis in Biometrics formulating IPMs as HMMs!
- Present CJS as a particular case of the multistate, as in https://github.com/oliviergimenez/multievent_jags_R
- Illustrate what dcat does; Categorical is a special case of the multinomial distribution with n = 1; The argument of the dcat distribution is a vector of probabilities for each category. 
- Prepare a slide or two on waic for cjs then anova, with sim and dipper
- Shall I use tablet/stylet?
- A few words to tell how we met back in Montpellier > 15 years ago Marc, Michael and I; Roger/Jean-Do in common, we're like cousin; introduce Maud 
- Tell how many neff we need for convergence

- ~~Add a few memes and animated gifs~~
- ~~Sort out what I will be demonstrating in the script~~
- ~~To maintain animations, duplicate slides instead of using powerpoint features, it makes the pdf conversion safer~~
- ~~Add bibliography w/ paper screenshots~~
- ~~Mention GOF tests w/ U-CARE and R2ucare, chapter on Gentle introduction~~
- ~~Add a slide w/ bibliography and/or To go further~~
- ~~Include pictures of key people in the capture-recapture world -> Rachel, Ruth, Eleni, Anita, etc.~~

On the practicals:

- Translate the R scripts in Rmd?
- ~~If too short to do the Rmd, just have a script for the lecture demo, and another one for exercises with and without solutions; show the matrices~~
- ~~Add an example of multievent model in the multistate R script? Use the roe deer example from my TPB paper https://github.com/oliviergimenez/multievent_jags_R; check out the examples we cover in our multievent workshops; make it an exercise;~~
- ~~The multistate example can be simplified and made faster, we don't need the loops on individuals/time to define the state matrices. In the exercise, we need the time loop. Mention that to add individual stuff, an extra loop on individuals is needed. Tension is with the lecture that presents the model with time and individual effect as the most general and generic model~~
- ~~Put m-array stuff in another file, gathering single and multistate stuff (or two files, yet to be decided)~~
- ~~Provide link to Nimble translation of book? Add Nimble implementation of the exercises? I won't do it; Maud has a full lecture to show-case Nimble~~
