



robu_Wild_parallel<- function(i){
  robu_Wild(formula = y_mat ~X5.10,
            data = Simdata_list[[i]],
            studynum = cluster,
            var.eff.size = Vars,
            userweights = weights,
            rho = 0.6,
            small = T,
            Replication_num = 1000,
            coef_index = 2)
}

Sys.time()
start_time <- Sys.time()
mc_result_single_X5.10 <- mclapply(X = 1:80, FUN = robu_Wild_parallel, mc.cores = detectCores())
end_time <- Sys.time()
time_interval <- end_time - start_time 

save(mc_result_single_X5.10, file="/Users/windshield/Desktop/Robumeta/Robu_Wild/Wild_boot_simulation/mc_result_single_X5.10.Rdata")



1 - sum(unlist(lapply(X = mc_result_single_X5.10,FUN = function(x){x[[3]]})) >0.05)/80

1 - sum(unlist(lapply(X = mc_result_single_X5.10,FUN = function(x){x[[1]]$reg_table[2,6]})) >0.05)/80






robu_Wild_parallel<- function(i){
  robu_Wild(formula = y_mat ~X5.10,
            data = Simdata_list[[i]],
            studynum = cluster,
            var.eff.size = Vars,
            userweights = weights,
            rho = 0.6,
            small = F,
            Replication_num = 1000,
            coef_index = 2)
}

Sys.time()
start_time <- Sys.time()
mc_result_single_X5.10_ns <- mclapply(X = 1:80, FUN = robu_Wild_parallel, mc.cores = detectCores())
end_time <- Sys.time()
time_interval <- end_time - start_time 

save(mc_result_single_X5.10_ns, file="/Users/windshield/Desktop/Robumeta/Robu_Wild/Wild_boot_simulation/mc_result_single_X5.10_ns.Rdata")



1 - sum(unlist(lapply(X = mc_result_single_X5.10_ns,FUN = function(x){x[[3]]})) >0.05)/80

1 - sum(unlist(lapply(X = mc_result_single_X5.10_ns,FUN = function(x){x[[1]]$reg_table[2,6]})) >0.05)/80






robu_Wild_parallel<- function(i){
  robu_Wild(formula = y_mat ~X1a,
            data = Simdata_list[[i]],
            studynum = cluster,
            var.eff.size = Vars,
            modelweights = "CORR",
            rho = 0.6,
            small = T,
            Replication_num = 1000,
            coef_index = 2)
}

Sys.time()
start_time <- Sys.time()
mc_result_single_X1a_corr <- mclapply(X = 1:80, FUN = robu_Wild_parallel, mc.cores = detectCores())
end_time <- Sys.time()
time_interval <- end_time - start_time 

save(mc_result_single_X1a_corr, file="/Users/windshield/Desktop/Robumeta/Robu_Wild/Wild_boot_simulation/mc_result_single_X1a_corr.Rdata")



1 - sum(unlist(lapply(X = mc_result_single_X1a_corr,FUN = function(x){x[[3]]})) >0.05)/80

1 - sum(unlist(lapply(X = mc_result_single_X1a_corr,FUN = function(x){x[[1]]$reg_table[2,6]})) >0.05)/80
