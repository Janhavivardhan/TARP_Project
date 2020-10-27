function [dn,sr,von1,voff1,vadtimes] = snreval(NOISY, varargin)

VERSION = 0.54;
DATE = 20140701;

% Parse out the optional arguments
[VAD,CLEAN,TS,TE,GUESS,DISP,LISTIN,LISTOUT, ...
 VADDIR,VADLIST,CLEANDIR,CLEANLIST,LDCLABELS,SAMPLERATE, ...
 CHECKFSHIFT, HPF, PREEMPH, MY_STNR, ...
 XTRA] = ...
    process_options(varargin, '-vad', '', '-clean', '', ...
                    '-start', 0, '-end', 0', '-guessvad', 0, ...
                    '-disp', 1, ...
                    '-listin', 0, '-listout', 0, ...
                    '-vaddir', '', '-vadlist', '', ...
                    '-cleandir', '', '-cleanlist', '', ...
                    '-ldclabels', 0, '-samplerate', 0, ...
                    '-checkfshift', 0, '-hpf', 0, '-preemph', 0, ...
                    '-my_stnr', 0);

HELP = 0;
if length(XTRA) > 0
  HELP = length(strmatch('-help',XTRA,'exact')) > 0;
  if ~HELP
    disp(['Unrecognized options:',sprintf(' %s',XTRA{1:end})]);
    HELP = 1;
  end
end

if strcmp(NOISY,'-help') == 1; HELP = 1; end
  
if HELP
  disp(['snreval v',num2str(VERSION),' of ',num2str(DATE)]);
  help('snreval');
  return
end

if LISTOUT || ~DISP
  QUIET = 1;
else
  QUIET = 0;
end

% 
% % Text header
% disp(['#============= SNREVAL v',num2str(VERSION),' (',num2str(DATE),') ===']);
% disp(['# args: ', join(varargin)]);
% % If LISTOUT, print out column headers
% if LISTOUT
%   disp(sprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s', ...
%                '#TARGETFILE', 't_start', 't_end', 't_delay', ...
%                'STNR', 'WADA', 'SNRvad', 'SAR', 'PESQ'));



do_guess_vad = GUESS;

% order of high-pass filters
HPFORD = 8;

SNRvad = 0;
SAD = 0;
pesqmos = 0;
targdelay = 0;

% default output, so doesn't crash if input list is empty
SNRstnr = [];

if LISTIN
  noisylist = listfileread(NOISY);
else
  noisylist = {NOISY};
end

if length(VADLIST) > 0
  vadlist = listfileread(VADLIST);
else
  vadlist = '';
end
if length(CLEANLIST) > 0
  cleanlist = listfileread(CLEANLIST);
else
  cleanlist = '';
end

nnoisy = length(noisylist);

for nf = 1:nnoisy
  
  % initialize values needed for list output
  targdelay = 0;
  SNRstnr = -999;
  SNRwada = -999;
  SNRvad = -999;
  SAR = -999;
  pesqmos = -999;

  
  NOISY = noisylist{nf};
  [p,n,e] = fileparts(NOISY);
  
  % maybe construct VAD and CLEAN files
  if length(VADDIR) > 0
    VAD = fullfile(VADDIR, [n,'.txt']);
  end
  if length(CLEANDIR) > 0
    CLEAN = fullfile(CLEANDIR, [n,e]);
  end
  if length(vadlist) > 0
    VAD = vadlist{nf};
  end
  if length(cleanlist) > 0
    CLEAN = cleanlist{nf};
  end

  % Read in the specified VAD file, if any
  if length(VAD) > 0
    [vadtimes,vnatimes] = read_vad_file(VAD, LDCLABELS);
    have_vad = 1;
  else
    have_vad = 0;
  end

  if length(CLEAN) > 0
    have_clean = 1;
  else
    have_clean = 0;
  end

  % Read in sound file
  DUR = max(0,TE-TS);
%   [dn,sr] = audioread(NOISY,SAMPLERATE,1,TS,DUR);
[dn,sr] = wavread(NOISY);
  if HPF > 0
    [bh,ah] = butter(HPFORD, HPF/(sr/2), 'high');
    dn = filter(bh, ah, dn);
  end
  if PREEMPH ~= 0
    dn = filter([1 -PREEMPH], 1, dn);
  end
  
  dlen = length(dn);

  [pp,nm,ee] = fileparts(NOISY);
  nm(nm=='_') = ' ';

  if DISP
    % Reset all the subplots
%     subplot(111)
%     plot(1)
  end
  
  if have_clean
%     [dc,src] = audioread(CLEAN,sr,1,TS,DUR);
[dc,src] = wavread(CLEAN);
    if HPF > 0
      [bh,ah] = butter(HPFORD, HPF/(sr/2), 'high');
      dc = filter(bh, ah, dc);
    end
    if PREEMPH ~= 0
      dc = filter([1 -PREEMPH], 1, dc);
    end
    dlen = min(dlen,length(dc));
  end

  if TE == 0 || TE > TS+dlen/sr
    thisTE = TS+dlen/sr;
  else
    thisTE = TE;
  end

  if thisTE > 0 || thisTE < TS+dlen/sr
    % truncate files as indicated
    dn = dn(1:round((thisTE-TS)*sr));
    dlen = length(dn);
    if have_clean
      dc = dc(1:round((thisTE-TS)*sr));
      dlen = min(dlen,length(dc));
    end
  end

  if have_vad == 0
    if do_guess_vad > 0
      if have_clean
        if ~QUIET
          disp(['Guessing VAD from clean file ', CLEAN,' ...']);
        end
        vadsource = dc;
      else
        if ~QUIET
          disp(['Guessing VAD from noisy file ', NOISY,' ...']);
        end
        vadsource = dn;
      end
      if do_guess_vad ~= 1.0
        vad_tsm = do_guess_vad;
      else
        vad_tsm = 0.25;
      end
      [von1,voff1] = guess_vad1(vadsource, sr, vad_tsm);
      vadtimes = TS + guess_vad(vadsource, sr, vad_tsm);
      vnatimes = [vadtimes(1:end-1,2),vadtimes(2:end,1)];
      have_vad = 1;
    else
      if ~QUIET
        disp('No VAD, but guessing is not selected');
      end
    end
  end
 
  % Choose fft size; 1024 for sr = 16000
  nfft = 2^round(log(1024 * sr/16000)/log(2));
  nhop = nfft/2;
  nsgframes = 1+ floor((dlen-nfft)/nhop);
  fr = sr/nhop;

  if have_vad
    % Convert VAD times to indices
    % sample level
    vv = zeros(1,length(dn))==1;
    % and also sgram frames (separate masks for active and inactive)
    vxf = zeros(1,nsgframes);
    nvf = zeros(1,nsgframes);
    lastvoff = 1+0;
    for i = 1:size(vadtimes,1)
      von = vadtimes(i,1);
      vof = vadtimes(i,2);
      % skip any entries entirely outside our time range
      if vof > TS && von < thisTE
        % map into our time range & clip
        von = max(0,von - TS);
        vof = min(thisTE-TS,vof - TS);
        vv((1+round(von*sr)):(round(vof*sr))) = (1==1);
        % only mark as voiced frames that are completely within voiced region
        vonf = von*fr;
        voff = vof*fr;
        vxf(ceil(1+vonf):floor(voff)) = 1;
      end
    end 
    for i = 1:size(vnatimes,1)
      non = vnatimes(i,1);
      nof = vnatimes(i,2);
      % skip any entries entirely outside our time range
      if nof > TS && non < thisTE
        % map into our time range & clip
        non = max(0,non - TS);
        nof = min(thisTE-TS,nof - TS);
        % only mark as unvoiced frames that are completely in unvoiced region
        nonf = non*fr;
        noff = nof*fr;
        nvf(ceil(1+nonf):floor(noff)) = 1;
      end
    end 
    % convert frame masks to indices
    vxf = find(vxf(1:nsgframes));
    nvf = find(nvf(1:nsgframes));
  else
    % no VAD
    vv = ones(1,length(dn))==1;  % to make a logical
    vxf = 1:nsgframes;
    nvf = [];
  end



  % Calculate NIST STNR
  % save the extracted region of the noisy file
%   SNRstnr = nist_stnr(dn,sr,MY_STNR);
%   msgstnr = ['NIST STNR = ',mynum2str(SNRstnr),' dB'];

  % Calculate the WADA SNR
%   SNRwada = wada_snr(dn,sr);
%   msgwada = ['WADA SNR  = ',mynum2str(SNRwada),' dB'];

  if DISP
    % Plot waveform and VAD
%     figure,
%     plot([1:length(dn)]/sr, dn);
    if have_vad
%       hold on; plot([0:length(vv)-1]/sr,0.2*vv,'-g'); hold off
%       hold on; plot([0:length(vv)-1]/sr,-0.2*vv,'-g'); hold off
    end
%     title(['File: ',nm,' - ',mynum2str(TS),'-',mynum2str(thisTE),' s']);
  end
end
end


function s = mynum2str(n)
s = sprintf('%.1f',n);
end

function s = join(array)
s = '';
for i = 1:length(array)
  val = array{i};
  if isstr(i)
    s = [s, val, ' '];
  else
    s = [s, num2str(val), ' '];
  end
end
end

