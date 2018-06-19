# 1. euclidian distance is not appropriate is muscles are independent. use unidimensional distance
# 2. # (MaxStep*(t>=10*dt) + 10.0*MaxStep*(t<10*dt)) = 2.0, which seems too big. 
	# You can't actually fit a vector of that size in the unit cube.

# Constraint values from a pdb trace session
# Brian Cohn June 12,2018 @bc
# (Pdb) U
# array([ 0.50704678,  0.64580454])
# (Pdb) Coefficient1, Coefficient2, Constraint1
# (379377442.8395462, -578605325.96008301, -179578150.00039524)
# (MaxStep*(t>=10*dt) + 10.0*MaxStep*(t<10*dt)) = 2.0

plot_2d_har_samples_with_prior_point <- function(samples_df, prior_point_2d_vec, bounds_tuple_of_numeric = list(c(0,1),c(0,1))){
	plot(samples_df[,1],
	 samples_df[,2],
	 xlim=bounds_tuple_of_numeric[[1]],
	 ylim=bounds_tuple_of_numeric[[2]], col="grey", pch=19, cex=0.08,
	 xlab="coefficient1", ylab="coefficient2")
	points(prior_point_2d_vec[1], prior_point_2d_vec[2], col="darkred", pch=21)
}

all_points_meet_equation <- function(samples, lefthand,righthand) {
	all(samples %>% 
	as.matrix %*% 
	lefthand - righthand < 1e-10
	)
}

#only for 2 dimensional systems
generate_constraint_via_AxByC_equation <- function(eqn_left_hand_side, eqn_right_hand_side)
{
	constr <- list(
	  constr = rbind(eqn_left_hand_side,-eqn_left_hand_side),
	  dir=rep('<=', 2),
	  rhs=c(eqn_right_hand_side,-eqn_right_hand_side))
	return(constr)
}

constraint_via_AxByC_with_bounds <- function(eqn_left_hand_side, eqn_right_hand_side, bounds_tuple_of_numeric){
	constr <- generate_constraint_via_AxByC_equation(eqn_left_hand_side, eqn_right_hand_side)
	constraint_object <- mergeConstraints(
	list(constr,
	# Variable 1
	upperBoundConstraint(2, 1, bounds_tuple_of_numeric[[1]][2]),
	lowerBoundConstraint(2, 1, bounds_tuple_of_numeric[[1]][1]),
	upperBoundConstraint(2, 2, bounds_tuple_of_numeric[[2]][2]),
	lowerBoundConstraint(2, 2, bounds_tuple_of_numeric[[2]][1]))
	)
return(constraint_object)
}

library(hitandrun)
# In the format of Ax + By = C, where A = coefficient1, B=coefficient2, C=constraint1, 
# and x and y represent the desired output variables
hit_and_run_delta_constrained <- function(max_delta_per_timestep,
											coefficient1,
											coefficient2,
											constraint1,
											n_samples,
											prior_point,
											bounds_tuple_of_numeric,
											plot=FALSE){

	eqn_left_hand_side <- c(1,coefficient1/coefficient2)
	eqn_right_hand_side <- constraint1/coefficient2

	constraint_object <- constraint_via_AxByC_with_bounds(eqn_left_hand_side, eqn_right_hand_side, bounds_tuple_of_numeric)
	state <- har.init(constraint_object)

	# Generate points without delta constraints
	samples <- har.run(state, n.samples=n_samples)$samples %>% as.data.frame
	colnames(samples) <- c("constraint1","constraint2")
	

	# Recompose with delta constraints
	constraint_object_with_delta_params <- mergeConstraints(list(
		constraint_object,
		upperBoundConstraint(2, 1, prior_point[1] + max_delta_per_timestep),
		lowerBoundConstraint(2, 1, prior_point[1] - max_delta_per_timestep),
		upperBoundConstraint(2, 2, prior_point[2] + max_delta_per_timestep),
		lowerBoundConstraint(2, 2, prior_point[2] - max_delta_per_timestep)
		)
	)
	state_delta <- har.init(constraint_object_with_delta_params, thin=100)
	samples_sliced <- har.run(state_delta, n.samples=n_samples)$samples %>% as.data.frame
	colnames(samples_sliced) <- c("constraint1","constraint2")
	
	# TODO write delta validator
	if (plot){
		plot_2d_har_samples_with_prior_point(samples, prior_point )
		points(samples_sliced, col="darkgreen", pch=19, cex=0.2)
	}
		valid <- all_points_meet_equation(samples_sliced,
		eqn_left_hand_side,
		eqn_right_hand_side)

	print(paste("All points from non-delta set are valid:", valid))
	return(samples_sliced)

}

hit_and_run_delta_constrained(max_delta_per_timestep = 0.2,
							  coefficient1 = 379377442.8395462,
							  coefficient2 = -578605325.96008301,
							  constraint1  = -179578150.00039524,
							  n_samples = 1000,
							  prior_point = c(0.50704678,0.64580454),
							  bounds_tuple_of_numeric = list(c(0,1),c(0,1)),
							  plot=FALSE)