function pars = get_sample_position(pars)

% ask for height and length of sample if not in Pars

if ~isfield(pars,'SampleL')

    pars.SampleL = input('Please give sample length in mm [default = 20]: ');
    if isempty(pars.SampleL)
        pars.SampleL = 20;
    end
end

if ~isfield(pars,'SampleH')

    pars.SampleH = input('Please give vertical position of sample center in the cavity,\nmeasured from the top collar in mm [default = 62.5]: ');
    if isempty(pars.SampleH)
        pars.SampleH = 62.5;
    end
end

end