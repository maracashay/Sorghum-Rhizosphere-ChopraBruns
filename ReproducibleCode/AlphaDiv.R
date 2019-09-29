##### 16. Alpha-Diversity Estimates and Plotting for 16S####

alpha.div<-estimate_richness(ps.rare, measures=c("Shannon", "Simpson", "Observed"))
alpha.div

shannon.model<-lmer(alpha.div$Shannon ~ Map.1$Line + (1|Map.1$Block))
anova(shannon.model) #(p = 0.1641)

simpson.model<-lmer(alpha.div$Simpson ~ Map.1$Line + (1|Map.1$Block))
anova(simpson.model) # p = .1916

observed.model<-lmer(alpha.div$Observed ~ Map.1$Line + (1|Map.1$Block))
anova(observed.model) # p = 0.1389


observed.model<-lmer(Observed ~ Frost 
                    * Line + (1|Block), data= alpha.div)
anova(observed.model)

plot(resid(observed.model), alpha.div$Observed)
qqmath(observed.model, id = 0.05)
# time results in a p < 0.05, indicating that the residuals are not normally distributed
# an run lme instead like so:
lme.obs <- lme(Observed~ Frost * Line, random = ~1|Block, data = alpha.div)
lme.obs.iden <- update(lme.obs,weights = varIdent(form=~1|Frost), data=alpha.div)

ms2<- dredge(lme.obs.iden)
subset(ms2, recalc.weights = FALSE)

# now let's do the same but for the model without varIdent
ms3<- dredge(lme.obs)
subset(ms3, recalc.weights=TRUE)

plot(lme.obs, resid(.) ~ fitted(.) | Frost)
anova(lme.obs, lme.obs.power)

residuals(lme.obs.power, type = "normalized")

# using AICc for small sample sizes
lme.full.varident<- lme(Observed~Frost*Line, random = ~1|Block, weights = varIdent(form=~1|Frost), data = alpha.div)

ms4<- dredge(lme.full.varident)
subset(ms4, recalc.weights = FALSE)

