import numpy as np

def bin2int(binArray):
  """
  Converts a numpy array to a binary number: bin2int(np.array([1,0,1,1])) = 13
  """
  tot = 0
  for j in range(len(binArray)):
    tot+=binArray[j]*(2**j)
  return int(tot)

def generateCodeMatrices(r,debug = True):
  if debug:
    print "Generating Hamming(%i,%i) code matrices"%(2**r-1,2**r-r-1)
  H = parityCheckMatrix(r)

  G = np.zeros((2**r-1,2**r-r-1))
  R = np.zeros((2**r-r-1,2**r-1))

  parityCols = [2**j-1 for j in range(r)]
  otherCols = [j for j in range(2**r-1) if j not in parityCols]

  parityRows = H[:,otherCols]
  G[parityCols,:] = parityRows
  G[otherCols,:] = np.eye(len(otherCols))

  R[:,otherCols] = np.eye(len(otherCols))
  return G,H,R

def parityCheckMatrix(r):
  H = []
  for j in range(1,(2**r)):
    H.append(map(int,list(bin(j))[2:])[-1::-1])
    H[-1].extend([0 for _ in range(r-len(H[-1]))])

  return np.array(H).T

class Message:
  """
  A class for playing with the Hamming codes.
  """
  def __init__(self,encodingBits = 3,mLength = 1000, debug = True):
    self.encodingBits = 3
    self.snippetLength = 2**self.encodingBits-self.encodingBits-1
    self.encodeSize = 2**self.encodingBits-1
    self.G,self.H,self.R = generateCodeMatrices(self.encodingBits,debug)
    self.generateMessage(mLength)
    self.encodedMessage = self.encodeMessage()
    self.debug = True

  def generateMessage(self,length):
    self.message = np.random.randint(0,2,self.snippetLength*(length/self.snippetLength))

  def encodeMessage(self):
    message = np.hstack((self.message,np.array([0 for _ in range(self.snippetLength-(len(self.message)-1)%self.snippetLength-1)])))
    message = message.reshape((message.size/self.snippetLength,self.snippetLength)).T
    encodedMessage = np.dot(self.G,message)
    return (encodedMessage.T).reshape(encodedMessage.size)%2

  def corruptMessage(self,nErrors):
    if self.debug:
      print 'Adding %i errors'%(nErrors)
    noise = np.zeros(self.encodedMessage.shape)
    indices = np.random.randint(0,self.encodedMessage.size,nErrors)
    noise[indices] = 1
    self.encodedMessage += noise

  def parityCheck(self):
    encodedMessage = self.encodedMessage.reshape((self.encodedMessage.size/self.encodeSize,self.encodeSize)).T
    parity = np.dot(self.H,encodedMessage)%2
    errorCols = parity.sum(0)
    if self.debug:
      print 'Found %i errors!'%(errorCols>0).sum()
    for col in range(len(errorCols)):
      if errorCols[col] > 0:
        row = bin2int(parity[:,col])-1
        encodedMessage[row,col] = (encodedMessage[row,col]+1)%2

  def decodeMessage(self):
    encodedMessage = self.encodedMessage.reshape((self.encodedMessage.size/self.encodeSize,self.encodeSize)).T
    decodedMessage = np.dot(self.R,encodedMessage)
    return (decodedMessage.T.reshape(decodedMessage.size))%2

def test(encodingBits = 3,mLength = 1000, nErrors = 10):
  M = Message(encodingBits,mLength)
  M.corruptMessage(nErrors)
  M.parityCheck()
  m = M.decodeMessage()
  if (m==M.message).all():
    print 'Decoded perfectly!'
  else:
    print 'Missed %i corrupted bits!'%((M.message!=m).sum())

if __name__ == '__main__':
  test()
