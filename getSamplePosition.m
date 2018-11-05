function Pars = getSamplePosition(Pars)

% ask for height and length of sample if not in Pars

if ~isfield(Pars,'SampleL') && ~isfield(Pars,'SampleH')

    Pars.SampleL = input('Please give sample length in mm [default = 20 mm]: ');
    if isempty(Pars.SampleL)
        Pars.SampleL = 20;
    end

    Pars.SampleH = input('Please give vertical position of sample center in cavity in mm,\nmeasured from the top collar [default = 62.5 mm]: ');
    if isempty(Pars.SampleH)
        Pars.SampleH = 62.5;
    end
end

end