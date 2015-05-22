function [ vlec_trellis ] = get_vlec_trellis( codewords, probabilites )
%GET_VLEC_TRELLIS Summary of this function goes here
%   Detailed explanation goes here

% codewords={[1 0 0 0 0 1]...
% [0 1 0 1 1 1]...
% [1 1 1 0 0]...
% [0 0 1 1 1]...
% [1 1 1 1]...
% [0 1 0 0]...
% [1 1 0]...
% [0 0 0]...
% [0 1 1]...
% [1 0 1]...
% [0 0 1 0]...
% [1 0 0 1]...
% [0 1 0 1 0]...
% [1 0 0 0 1]...
% [0 0 1 1 0 0]...
% [1 1 1 0 1 0]};
% 
% probabilites = ones(1,numel(codewords))/numel(codewords);

max_len = 0;
for i = 1:numel(codewords)   
    if numel(cell2mat(codewords(i))) > max_len;
        max_len = numel(cell2mat(codewords(i)));
    end
    
end

num_codewords = numel(codewords);

trans_from = ones(num_codewords,max_len)*-1;
trans_to = ones(num_codewords,max_len)*-1;
cw_table = ones(num_codewords,max_len)*-1;

trans_codeword_map = zeros(1,num_codewords*max_len);  %max length, will be mostly unused
state_probs = zeros(1,num_codewords*max_len);  %max length, will be mostly unused

next_state = 1;
tran_count = 0;
% % %initialise with the first entry
% % c = cell2mat(codewords(1));
% % for j = 1:length(c)
% %     cw_table(1,j) = c(j);
% %     trans_from(1,j) = next_state;
% %     if j == length(c)  %last transition goes back to 1
% %          trans_to(1,j) = 1;
% %     else                
% %         next_state = next_state + 1;
% %         trans_to(1,j) = next_state;
% %         tran_count = tran_count + 1;
% %         trans_codeword_map(tran_count) = i;
% %     end
% %     
% %     state_probs( trans_from(1,j) ) = state_probs( trans_from(1,j) ) + probabilites(1);
% % end
% % %state_probs( trans_to(1,j) ) = state_probs( trans_to(1,j) ) + probabilites(1);

%now do the rest
for i = 1:numel(codewords)
   
    %get our codeword in an array
    c = cell2mat(codewords(i));    
    
    %look for which existing entry has most similar prefix
    prelen = 0;
    preindex = 0;
    
    for k = 1:i
        pl = prefix_len(c,cw_table(k,:));
        if (pl > prelen)
            prelen = pl;
            preindex = k;            
        end
    end
    
    if (prelen > 0)
        state_from =  trans_to(preindex,prelen);
    else
        state_from = 1;
    end
    
    if length(c) <= prelen
       error('Not valid codewords'); 
    end
    
    for j = 1:length(c)        
        
        cw_table(i,j) = c(j);
        
        if j > prelen  %new transition
            trans_from(i,j) = state_from;
            tran_count = tran_count + 1;
            
            if j == length(c)  %last transition goes back to 1
               trans_to(i,j) = 1;
               trans_codeword_map(tran_count) = i;
            else
                next_state = next_state + 1;
                trans_to(i,j) = next_state;
                state_from = next_state;
                
            end            
            
            
            
            state_probs( trans_from(i,j) ) = state_probs( trans_from(i,j) ) + probabilites(i);
        else  %existing transition
            trans_to(i,j) = trans_to(preindex,j);
            trans_from(i,j) = trans_from(preindex,j);
            state_probs( trans_from(preindex,j) ) = state_probs( trans_from(preindex,j) ) + probabilites(i);
        end
    end
    %state_probs( trans_to(i,j) ) = state_probs( trans_to(i,j) ) + probabilites(i);
    
end

%consolidate everything into one matrix

vlec_trellis = zeros(tran_count,5);
%state from, state to, uncoded, coded,  log prob

tran = 1;
for i = 1:numel(codewords)
    for j = 1:max_len
        
        if trans_to(i,j) > 0
            fr = trans_from(i,j);
            to = trans_to(i,j);
            un = cw_table(i,j);
            %en = ????????????
            
            if findpair(fr,to,vlec_trellis(:,1:2)) == 0  %if a new transition
                vlec_trellis(tran,1) = fr;
                vlec_trellis(tran,2) = to;
                vlec_trellis(tran,3) = un;
                if (to > 1)
                    vlec_trellis(tran,5) = log(state_probs(to)/state_probs(fr));
                else
                    c = trans_codeword_map(tran);
                    vlec_trellis(tran,5) = log( probabilites(c)/state_probs(fr));
                end
                tran = tran + 1;
            end
            
        end
        
    end
    
end


end

function [found] = findpair(val1, val2, matrix_in)
found = 0;
for i = 1:size(matrix_in,1);
   if ((matrix_in(i,1) == val1) &&  (matrix_in(i,2) == val2))
      found = 1;
      return;
   end
end

end

function [prefix_len] = prefix_len(input1, input2)
    len = max(numel(input1),numel(input2));
    prefix_len = 0;
    
    for i = 1:len
       if input1(i) ~= input2(i)
           return;
       else
           prefix_len = prefix_len + 1;
       end
    end
        
end