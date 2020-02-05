function coef = R2(y_raw, y_fit)

coef=1-sum((y_raw - y_fit).^2)/sum((y_raw-mean(y_raw)).^2);

end