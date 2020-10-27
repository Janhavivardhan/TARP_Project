filename = 'user1.wav' % Calculate the codebook vector.
    for k=1:5
    for j=1:7
    cep = learn(['data\',num2str(k),'\',num2str(j),'.wav']);
     temp= learn(filename);
     if size(cep,1)>size(temp,1)
     cep=cep(1:size(temp,1),:);
     else
         temp=temp(1:size(cep,1),:);
     end
%     codebook = vqlbg(cep, 8);
    result(j) = corr(cep(:),temp(:));%identify(filename, codebook);
    end
    out(k)=max(result);
    end
    out