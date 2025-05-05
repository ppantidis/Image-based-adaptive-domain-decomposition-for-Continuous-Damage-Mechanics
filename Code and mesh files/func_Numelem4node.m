function[numelem4node]=func_Numelem4node(nnodes, connect)

numelem4node = zeros(nnodes,1);

for i=1:size(connect,1)
    for j=1:size(connect,2)
        numelem4node(connect(i,j))=numelem4node(connect(i,j))+1;
    end
end

end