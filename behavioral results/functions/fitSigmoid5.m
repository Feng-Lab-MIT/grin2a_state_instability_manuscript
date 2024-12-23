function Y_model = fitSigmoid5(BinaryData)
    type = 'Normal';
    model = @(k,x)  (1-cdf(type,x,k(1),k(2)));
    k = nlinfit([1:length(BinaryData)]',smoothdata(BinaryData,'gaussian',6,'omitnan'),model,[length(BinaryData)/2 1],'Options',statset('FunValCheck','off','MaxIter',20000,'WgtFun','logistic'));
    k(2) = k(2)/3;
    Y_model = model(k,[1:length(BinaryData)]');
end