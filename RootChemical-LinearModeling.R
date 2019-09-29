X3da.test<-lmer(X3DAs ~ Frost * Line + (1|Block), map.1)
anova(X3da.test) 

leveneTest(residuals(X3da.test) ~ map.1$Frost)
##### Stat Difference in DPPH ####
dpph.test<-lmer(DPPH ~ Frost * Line + (1|Block), map.1)
anova(dpph.test) # interaction was significant

leveneTest(residuals(dpph.test) ~ map.1$Frost)
leveneTest(residuals(dpph.test)~ map.1$Line)

em.dpph <-emmeans(dpph.test, ~ Line)
pairs.em <- pairs(em.dpph, simple = "each", combine = TRUE, p.adjust.methods= "fdr")
pairs.em

