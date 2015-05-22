addpath('testing');

inc = randn([1 1000]);
    inu = randn([1 1000]);
[c,d]=URC8_decoder_dual(inc,inu);
[a,b]=bcjr_decoder(inu,inc,'urc8');

if (sum(abs(a-c)) + sum(abs(b-d))>0.1)
    error('somethings wrong 1');
end


transitions_uec = gen_uec_transitions( [1;1;1;1], [0.797000539514745,0.116806964340725,0.0379853325346553], 1.5368 );
[g,h]  = bcjr_decoder(zeros(1,1000),inc,'uec',transitions_uec);
[t, y] = trellis_decoder(inc, [1;1;1;1], [0.797000539514745,0.116806964340725,0.0379853325346553], 1.5368);

inc2 = randn([2 1000]);
transitions_uec1 = gen_uec_transitions([0 1; 1 1], 0.797000539514745, 1.5368 );
[g1,h1]  = bcjr_decoder(zeros(1,1000),inc2,'uec',transitions_uec1);
[t1, y1] = trellis_decoder(inc2, [0 1; 1 1], 0.797000539514745, 1.5368);


inc2 = randn([2 1000]);
trans2 =  [1,          1,          0,          0,          0,  0;
                1,          3,          1,          1,          1,  0;
                2,          1,          0,          0,          1,  0;
                2,          3,          1,          1,          0,  0;
                3,          4,          0,          1,          0,  0;
                3,          2,          1,          0,          1,  0;
                4,          4,          0,          1,          1,  0;
                4,          2,          1,          0,          0,  0];
[un21,en21]  = bcjr_decoder(zeros(1,1000),inc2,'uec',trans2);
[un22,en22] = CC2_decoder_bcjr(zeros(1,1000), inc2);


inc3 = randn([3 1000]);
trans3 =       [1,          1,          0,          0,          0,           0,  0;
                1,          3,          1,          1,          1,           1,  0;
                2,          1,          0,          0,          1,           1,  0;
                2,          3,          1,          1,          0,           0,  0;
                3,          4,          0,          1,          0,           0,  0;
                3,          2,          1,          0,          1,           1,  0;
                4,          4,          0,          1,          1,           1,  0;
                4,          2,          1,          0,          0,           0,  0];
[un31,en31]  = bcjr_decoder(zeros(1,1000),inc3,'uec',trans3);
[un32,en32] = CC3_decoder_bcjr(zeros(1,1000), inc3);


feval('zeta_dist_p1_07_N_10_211')
vlec_codebook = VLCTable;
vlec_probs = SymbolProbs;
vlec_trellis = get_vlec_trellis( vlec_codebook, vlec_probs );
w_a = randn([1 1000]);

unv1 = vlec_viterbi_decoder(w_a, vlec_trellis);
unv2 = viterbi_decoder(w_a, vlec_trellis);


if (sum(abs(t-h)) + sum(abs(y-g))>0.1)
    error('somethings wrong 2');
end
if (sum(sum(abs(t1-h1))) + sum(abs(y1-g1))>0.1)
    error('somethings wrong 3');
end
if (sum(sum(abs(un21-un22))) + sum(abs(en21-en22))>0.1)
    error('somethings wrong 4');
end
if (sum(sum(abs(un31-un32))) + sum(abs(en31-en32))>0.1)
    error('somethings wrong 5');
end
if (sum(sum(abs(unv1-unv2)))>0.1)
    error('somethings wrong 6');
end


%%%%%%% URC8

tic

for i=1:100
    inc = randn([1 1000]);
    inu = randn([1 1000]);

    [c,d]=URC8_decoder_dual(inc,inu);
    
end

toc



tic

for i=1:100
    inc = randn([1 1000]);
    inu = randn([1 1000]);

    [c,d]=bcjr_decoder(inu,inc,'urc8');
    
end

toc

%%%%%%% UEC 1 codewords
tic

for i=1:100
    inc = randn([1 1000]);
    inu = randn([1 1000]);

    [c,d]=trellis_decoder(inc, [1;1;1;1], [0.797000539514745,0.116806964340725,0.0379853325346553], 1.5368);
    
end

toc



tic

for i=1:100
    inc = randn([1 1000]);
    inu = randn([1 1000]);

    [c,d]=bcjr_decoder(zeros(1,1000),inc,'uec',transitions_uec);
    
end

toc

%%%%%%% UEC 2 codewords
tic

for i=1:100
    inc = randn([2 1000]);
    inu = randn([1 1000]);

    [c,d]=trellis_decoder(inc, [0 1; 1 1], [0.797000539514745], 1.5368);
    
end

toc



tic

for i=1:100
    inc = randn([2 1000]);
    inu = randn([1 1000]);

    [c,d]=bcjr_decoder(zeros(1,1000),inc,'uec',transitions_uec1);
    
end

toc

%%%%%%% VLEC codewords

tic

for i=1:100
    inu = randn([1 1000]);
    unv1 = vlec_viterbi_decoder(inu, vlec_trellis);
    
end

toc
tic

for i=1:100

    inu = randn([1 1000]);    
    unv2 = viterbi_decoder(inu, vlec_trellis);
    
end

toc
