#List to save all warning messages
warnings_list <- vector('list',2)


#First fit with a single warning
dat <- data.frame(
    x = c(1, 2, 3, 4, 5, 6),
    y = c(0, 0, 0, 1, 1, 1)
)


fit <- withCallingHandlers(
    glm(y ~ x, family = binomial(), data = dat),
    warning = function(w) {
        warnings_list[[1]] <<- c(warnings_list[[1]] , conditionMessage(w))
    }
)




#Second fit with two warnings (note that I have purposefully set maxit=10 to create two warning messages)
dat2 <- data.frame(
    x = c(1, 2, 3, 4, 5, 6),
    y = c(0, 0, 0, 1, 1, 1)
)


fit2 <- withCallingHandlers(
    glm(y ~ x, family = binomial(), data = dat2, control = glm.control(maxit = 10)),
    warning = function(w) {
      warnings_list[[2]] <<- c(warnings_list[[2]] , conditionMessage(w))
    }
)


#End result is a list containing warning message(s) for each fit
warnings_list
