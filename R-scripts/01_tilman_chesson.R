
library(tidyverse)
library(cowplot)

## get data
all_rstars <- read_csv("data-processed/all-rstars.csv")
all_combos <- read_csv("data-processed/all_combos_allopatry.csv") %>% 
	split (.$combination)

### set up Tilman to Chesson conversion, including mapping from supply point to zones for alpha conversion
## first set up some supply points
ns <- c(seq(1*(min(all_rstars$n_star) - 0.01), 1*(max(all_rstars$n_star) + 0.01), length.out = 4), 
		seq(2*(min(all_rstars$n_star) - 0.01), 2*(max(all_rstars$n_star) + 0.01), length.out = 4),
		seq(3*(min(all_rstars$n_star) - 0.01), 3*(max(all_rstars$n_star) + 0.01), length.out = 4))
ps <- c(seq(1*(min(all_rstars$p_star) - 0.001), 1*(max(all_rstars$p_star) + 0.001), length.out = 4),
		seq(2*(min(all_rstars$p_star) - 0.001), 2*(max(all_rstars$p_star) + 0.001), length.out = 4),
		seq(3*(min(all_rstars$p_star) - 0.001), 3*(max(all_rstars$p_star) + 0.001), length.out = 4))


### Tilman to Chesson mechanics

find_alphas_combos <- function(combos, SN1, SP1){
	snippet <- all_rstars %>% 
		filter(population %in% c(combos$comp1[[1]], combos$comp2[[1]])) 
	
	pop1 <- snippet$population[[1]]
	pop2 <- snippet$population[[2]]
	
	c1P <- snippet$pc[snippet$population == pop1]
	c2P <- snippet$pc[snippet$population == pop2]
	c1N <- snippet$nc[snippet$population == pop1]
	c2N <- snippet$nc[snippet$population == pop2]
	R1P <- snippet$p_star[snippet$population == pop1]
	R2N <- snippet$n_star[snippet$population == pop2]
	R2P <- snippet$p_star[snippet$population == pop2]
	R1N <- snippet$n_star[snippet$population == pop1]
	D <- 0.5 ## set dilution
	
	## define the consumption vector lines
	
	SN <- SN1
	SP <- SP1
	
	
	r1 <- max(snippet$p_umax[snippet$population == pop1], snippet$n_umax[snippet$population == pop1])
	r2 <- max(snippet$p_umax[snippet$population == pop2], snippet$n_umax[snippet$population == pop2])
	
	cons_vec1_intercept <- max(R1P, R2P) -(c1P/c1N)*max(R1N, R2N)
	cons_vec2_intercept <- max(R1P, R2P) -(c2P/c2N)*max(R1N, R2N)
	supply_vec <- SP/SN
	
	cons_vec1_fun <- function(x){
		y <- (c1P/c1N)*x + cons_vec1_intercept
		return(y)
	}
	
	cons_vec2_fun <- function(x){
		y <- (c2P/c2N)*x + cons_vec2_intercept
		return(y)
	}
	
	supply_vec_fun <- function(x){
		y <- (SP/SN)*x
		return(y)
	}
	
	zone_middle <- cons_vec2_fun(SN) <= supply_vec_fun(SN) & supply_vec_fun(SN) <= cons_vec1_fun(SN) | cons_vec1_fun(SN) <= supply_vec_fun(SN) & supply_vec_fun(SN) <= cons_vec2_fun(SN)
	zone_bottom <- cons_vec2_fun(SN) >= supply_vec_fun(SN) & supply_vec_fun(SN) <= cons_vec1_fun(SN)
	zone_top <- cons_vec2_fun(SN) <= supply_vec_fun(SN) & supply_vec_fun(SN) >= cons_vec1_fun(SN)
	
	zones <- c("zone_middle", "zone_bottom", "zone_top")
	zone <- zones[c(zone_middle, zone_bottom, zone_top)]
	
	alphas <- function(zone) {
		if (zone == "zone_middle") {
			a11 <- c1P / (D * (SP - R1P))
			a12 <- c2P / (D * (SP - R1P))
			a21 <- c1N / (D * (SN - R2N))
			a22 <- c2N / (D * (SN - R2N))
		} else if (zone == "zone_top") {
			a11 <- c1N / (D * (SN - R1N))
			a12 <- c2N / (D * (SN - R1N))
			a21 <- c1N / (D * (SN - R2N))
			a22 <- c2N / (D * (SN - R2N))
		} else if (zone == "zone_bottom") {
			a11 <- c1P / (D * (SP - R1P))
			a12 <- c2P / (D * (SP - R1P))
			a21 <- c1P / (D * (SP - R2P))
			a22 <- c2P / (D * (SP - R2P))
		}
		alphas1 <- data.frame(a11 = a11, a12 = a12, a21 = a21, a22 = a22)
		return(alphas1)
	}
	
	alphas_calc <- alphas(zone)
	r_s <- data.frame(r1 = r1, r2 = r2)
	comps <- data.frame(pop1 = pop1, pop2 = pop2, D = D, treatment = snippet$treatment[1], nitrogen_supply = SN1,
						phosphorus_supply = SP1, zone = zone)
	alphas_calc2 <- bind_cols(alphas_calc, r_s, comps)
	return(alphas_calc2)
}


## loop through some resource supply points that span the resource space from below to 3x above the R*s
## ns and ps are the different amounts of N and P supply
results <- data.frame()
for(i in ns){
	for(j in ps){
		hold <- map_df(all_combos, find_alphas_combos, SN1 = i, SP1 = j)
		hold$n_supply <- i
		hold$p_supply <- j
		results <- bind_rows(results, hold)
	}}

results2b <- results %>% 
	mutate(rho = sqrt((a12*a21)/(a11*a22))) %>% 
	mutate(fit_ratio = (sqrt((a11*a12)/(a22*a21)) * (r2-D)/(r1-D))) %>% 
	mutate(coexist = rho < fit_ratio &  fit_ratio < 1/rho) %>% 
	mutate(stabil_potential = 1 - rho) %>% 
	mutate(supply_ratio = n_supply/p_supply)

### ok this looks right, for the combination of population 10 and 1, and our three supply points, with one chosen supply ratio
### zngis are crossing and zones are mapping correctly

snip <- results2b %>%
	mutate(combination = paste(pop1, pop2, sep = "_")) %>% 
	group_by(supply_ratio) %>% 
	mutate(supply_ratio_round = round(supply_ratio, digits = 3)) %>% 
	filter(combination == "10_1") %>% 
	filter(supply_ratio_round == 10.681) 

all_rstars %>% 
	filter(population %in% c(10, 1)) %>% 
	ggplot(aes(x = n_star, y = p_star)) + geom_point(color = "grey") +
	geom_segment(aes(x = max(n_star), y = max(p_star), xend = n_star + nc*9, yend = p_star + pc*9),
				 size = 0.5, linetype = 2, color = "grey") +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star, yend = 1),
				 size = 0.5, linetype = 1, color = "grey") +
	geom_segment(aes(x = n_star, y = p_star, xend = 4, yend = p_star),
				 size = 0.5, linetype = 1, color = "grey") +
	ylab("P (uM)") + xlab("N (uM)") +
	geom_point(aes(x = nitrogen_supply, y = phosphorus_supply, color = zone), data = snip) +
	coord_cartesian() +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black") 
ggsave("figures/population_10_1_combo.png", width = 8, height = 6)

## this is clearly wrong, zngis aren't crossing, and the consumption vectors don't cross at the right place
## because the consumption vectors are written such they cross at the elbows of the crossing zngis


snip2 <- results2b %>%
	mutate(combination = paste(pop1, pop2, sep = "_")) %>% 
	group_by(supply_ratio) %>% 
	mutate(supply_ratio_round = round(supply_ratio, digits = 3)) %>% 
	filter(combination == "22_6") %>% 
	filter(supply_ratio_round == 7.325) 

all_rstars %>% 
	filter(population %in% c(22, 6)) %>% 
	ggplot(aes(x = n_star, y = p_star)) + geom_point(color = "grey") +
	geom_segment(aes(x = max(n_star), y = max(p_star), xend = n_star + nc*9, yend = p_star + pc*9),
				 size = 0.5, linetype = 2, color = "grey") +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star, yend = 1),
				 size = 0.5, linetype = 1, color = "grey") +
	geom_segment(aes(x = n_star, y = p_star, xend = 4, yend = p_star),
				 size = 0.5, linetype = 1, color = "grey") +
	ylab("P (uM)") + xlab("N (uM)") +
	geom_point(aes(x = nitrogen_supply, y = phosphorus_supply, color = zone), data = snip) +
	coord_cartesian() +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black") 
ggsave("figures/population_22_6_combo.png", width = 8, height = 6)

