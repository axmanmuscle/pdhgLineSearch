
function [out, phaseImg] = partial_fourier_b( in, sFSR, varargin )
  % out = mri_reconPartialFourier( in, sFSR [, 'phases', phases, 'op', op ] )
  %
  % Written according to "Partial k-space Reconstruction" by John Pauly
  %
  % Inputs:
  % in - the input array of size Ny x Nx x nCoils representing the MRI data.
  % sFSR - Either a scalar or a two element array that specifies the size of the fully
  %   sampled region, which is used to estimate the phase of the image.
  %   If a scalar, then the size of the fully sampled region is sFSR x Nx.
  %   If a two element array, then the size of the FSR is sFSR(1) x sFSR(2).
  %
  % Outputs:
  % out - an array of size Ny x Nx x nCoils that represents the coil images
  %
  % Written by Nicholas Dwork - Copyright 2023
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'op', 'notransp', @(x) true );
  p.addParameter( 'phases', [], @isnumeric );
  p.parse( varargin{:} );
  phases = p.Results.phases;
  op = p.Results.op;

  Ny = size( in, 1 );
  ys = size2imgCoordinates( Ny );
  centerIndx = find( ys == 0, 1 );

  if numel( sFSR ) == 1
    sFSR = [ sFSR, size( in, 2 ) ];
  end

  firstY = centerIndx - ceil( ( sFSR(1) - 1 ) / 2 );
  lastY = centerIndx + floor( ( sFSR(1) - 1 ) / 2 );
  m = 2 / ( lastY - firstY );
  ramp = m * ys + 1.0;
  ramp( 1 : firstY ) = 0;
  ramp( lastY : end ) = 2;

  nCoils = size( in, 3 );

  if numel( phases ) == 0
    if strcmp( op, 'transp' )
      error( 'Cannot estimate image phase during transpose operation' );
    end

    firstX = centerIndx - ceil( ( sFSR(2) - 1 ) / 2 );
    lastX = centerIndx + floor( ( sFSR(2) - 1 ) / 2 );

    inLF = in;  % Low-freq data
    if nCoils == 1
      inLF( 1:firstY-1, : ) = 0;
      inLF( :, 1:firstX-1 ) = 0;
      inLF( :, lastX+1:end ) = 0;
    else
      inLF( 1:firstY-1, :, : ) = 0;
      inLF( :, 1:firstX-1, : ) = 0;
      inLF( :, lastX+1:end, : ) = 0;
    end
    imgLF = fftshift2( ifft2( ifftshift2( inLF ) ) );
    phases = angle( imgLF );
  
  else
    phaseImg = exp( 1i * phases );
    %%% Form the explicit operator for B and see that you get the desired
    %%% outcome
    %%% AA^T + BB^T = 1/alpha I
  
    %%% For testing just for Ax
    inW = bsxfun( @times, in, 2-ramp );
    imgW = fftshift2( ifft2( ifftshift2( inW ) ) );
    tmp1 = imgW .* conj( phaseImg );
    Ax = real( tmp1 ) .* phaseImg;

    %%% Form A^Tx
    Atx = real( in .* conj( phaseImg ) );
    inWithPhase = Atx .* phaseImg;
    ifft2hIn = fftshift2( ifft2h( ifftshift2( inWithPhase ) ) );
    Atx = bsxfun( @times, ifft2hIn, 2-ramp );

    inW2 = bsxfun( @times, Atx, 2-ramp );
    imgW2 = fftshift2( ifft2( ifftshift2( inW2 ) ) );
    tmp12 = imgW2 .* conj( phaseImg );
    AAtx = real( tmp12 ) .* phaseImg;


    %%% split apart and apply the ramp
    n = size(in, 1);
    insplit = [real(in) ; imag(in)];
    tmpramp = 2 - ramp;
    splitramp = [diag(tmpramp) zeros(n); zeros(n) diag(tmpramp)];
    tmp2 = splitramp * insplit;
    test2 = tmp2(1:n, :) + 1i*tmp2(n+1:end, :);

    norm(test2 - inW);

    %%% Try doing AA^Tx with the split ramp^2
    isplit = [real(ifft2hIn); imag(ifft2hIn)];
    tmp3 = splitramp.^2 * isplit;
    out3 = tmp3(1:n, :) + 1i*tmp3(n+1:end, :);
    out4 = fftshift2( ifft2( ifftshift2( out3 ) ) );
    tmp4 = out4 .* conj( phaseImg );
    AAtx2 = real( tmp4 ) .* phaseImg;

    norm(AAtx - AAtx2)


    %%% form B^Tx
    alpha = 1/8;
    n2 = zeros(size(splitramp));
    Bramp = sqrt(1/alpha * n2 - splitramp.^2);

    Btx = real( in .* conj( phaseImg ) );
    inWithPhase = Btx .* phaseImg;
    ifft2hIn = fftshift2( ifft2h( ifftshift2( inWithPhase ) ) );
    isplit = [real(ifft2hIn); imag(ifft2hIn)];
    out5 = Bramp * isplit;
    out6 = out5(1:n, :) + 1i*out5(n+1:end, :);
    out7 = fftshift2( ifft2( ifftshift2( out6 ) ) );
    tmp7 = out7 .* conj( phaseImg );
    BBtx2 = real( tmp7 ) .* phaseImg;

    norm(AAtx2 + BBtx2 - x/alpha)
    


    

  end
  phaseImg = exp( 1i * phases );

  if strcmp( op, 'notransp' )
    inW = bsxfun( @times, in, 2-ramp );
    imgW = fftshift2( ifft2( ifftshift2( inW ) ) );
    out = imgW .* conj( phaseImg );
    out = real( out ) .* phaseImg;

  elseif strcmp( op, 'transp' )
    in = real( in .* conj( phaseImg ) );
    inWithPhase = in .* phaseImg;
    ifft2hIn = fftshift2( ifft2h( ifftshift2( inWithPhase ) ) );
    out = bsxfun( @times, ifft2hIn, 2-ramp );

  else
    error( 'Unrecognized operator' );
  end
end
