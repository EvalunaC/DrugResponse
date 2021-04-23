rf.model <- randomForest(Resp~.,data=trainFrame,
                         ntree = 50,
                         nodesize = 5,
                         mtry = 2,
                         importance = TRUE,
                         metric = "RMSE")

caret.oob.model <- train(Resp~.,data=trainFrame,
                         method = "rf",
                         ntree = 50,
                         tuneGrid = data.frame(mtry = 2),
                         nodesize = 5,
                         importance = TRUE,
                         metric = "RMSE",
                         trControl = trainControl(method = "oob", seed = c(1,1)),
                         allowParallel = FALSE)

seeds <- as.vector(c(1:26), mode = "list")

# For the final model
seeds[[26]] <- 1

caret.boot.model <- train(Resp~.,data=trainFrame,
                          method = "rf",
                          ntree = 50,
                          tuneGrid = data.frame(mtry = 2),
                          nodesize = 5,
                          importance = TRUE,
                          metric = "RMSE",
                          trControl = trainControl(method = "boot", seeds = seeds),
                          allowParallel = FALSE)

rf.model
caret.oob.model$final
caret.boot.model$final

best <- as.numeric(row.names(caret.oob.model$bestTune))
rf_result <- list(method = "Random Forest",
                   RMSE = caret.oob.model$result$RMSE[best],
                   R_Square = caret.oob.model$result$Rsquared[best],
                   MAE = caret.oob.model$results$MAE[best]
                  )
caret_rf <- predict(caret.oob.model,trainFrame)
eval_result(caret_rf)
postResample(pred = caret_rf, obs = trainFrame$Resp)

eval_result<-function(preds){
  MAE <- mean(abs(trainFrame$Resp - preds))
  SSE <- sum((trainFrame$Resp - preds)^2)
  SST <- sum((trainFrame$Resp - mean(trainFrame$Resp))^2)
  SSM <- sum((preds-mean(trainFrame$Resp))^2)
  R_square <- 1 - SSE / SST
  R2 <- (cor(preds,trainFrame$Resp))^2
  R2adj<-1-((1-R_square)*(nrow(trainFrame)-1)/(nrow(trainFrame)-(ncol(trainFrame)-1)-1))
  RMSE = sqrt(SSE/nrow(trainFrame))
  F_stat<-SSM/(ncol(trainFrame)-1)/(SSE/(nrow(trainFrame)-ncol(trainFrame)))
  t_test<-t.test(trainFrame$Resp, preds)$p.value
  ks_test<-ks.test(trainFrame$Resp, preds)$p.value
  results<-list(R2=R2,RMSE=RMSE,R_Square=R_square, Adjusted_R2=R2adj,MAE=MAE, F_stat=F_stat,t_test=t_test,ks_test=ks_test)
  return(results)
}
